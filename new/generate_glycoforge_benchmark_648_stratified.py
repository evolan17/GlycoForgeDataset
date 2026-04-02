from __future__ import annotations

import argparse
import csv
import inspect
import json
import random
import shutil
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.tools.registry import export_tool_specs_for_manifest  # noqa: E402

DEFAULT_OUT_ROOT = Path("data/glycoforge_benchmark_648")
DEFAULT_FALLBACK_SOURCE_PACK = Path("data/glycoforge_benchmark_500")
BASE_SEED = 20260401

TASK_FAMILIES = ["single_workflow", "batch_correction_centric", "paired_coupled"]
COMPLEXITY_TIERS = ["simple", "moderate", "hard"]
SIGNAL_REGIMES = ["absent_or_weak", "moderate", "strong_or_confounded"]
MODE_LEVELS = ["synthetic", "templated"]
MISSING_LEVELS = ["low", "medium", "high"]
BATCH_LEVELS = ["light", "elevated"]
MOTIF_LEVELS = [False, True]
PER_CLUSTER_DEFAULT = 24


@dataclass(frozen=True)
class CaseBlueprint:
    base_mode: str
    missing_level: str
    batch_level: str
    motif_driven: bool


@dataclass(frozen=True)
class CaseSpec:
    dataset_idx: int
    dataset_id: str
    cluster_id: str
    task_family: str
    complexity_tier: str
    signal_regime: str
    blueprint_index: int
    base_mode: str
    missing_level: str
    batch_level: str
    motif_driven: bool
    paired_requested: bool
    requested_mode: str
    seed: int


@dataclass
class SourceCase:
    dataset_id: str
    dataset_dir: Path
    manifest: dict[str, Any]
    hidden_truth: dict[str, Any]
    task_card: dict[str, Any] | None
    missing_fraction: float
    batch_strength: float
    bio_strength: float
    source_mode: str
    source_scenario: str
    has_missing: bool
    has_batch: bool


def _try_import_simulate():
    try:
        from glycoforge import simulate
    except Exception:
        return None
    return simulate


def _try_import_simulate_paired():
    try:
        from glycoforge import simulate_paired
    except Exception:
        return None
    return simulate_paired


def _simulate_param_names(simulate_fn) -> set[str]:
    try:
        return set(inspect.signature(simulate_fn).parameters)
    except Exception:
        return set()


def _build_blueprints() -> list[CaseBlueprint]:
    blueprints = [
        CaseBlueprint(base_mode=mode, missing_level=missing, batch_level=batch, motif_driven=motif)
        for mode, missing, batch, motif in product(MODE_LEVELS, MISSING_LEVELS, BATCH_LEVELS, MOTIF_LEVELS)
    ]
    blueprints.sort(
        key=lambda bp: (
            MODE_LEVELS.index(bp.base_mode),
            MISSING_LEVELS.index(bp.missing_level),
            BATCH_LEVELS.index(bp.batch_level),
            int(bp.motif_driven),
        )
    )
    return blueprints


def _build_case_specs(per_cluster: int) -> list[CaseSpec]:
    if per_cluster < 1 or per_cluster > PER_CLUSTER_DEFAULT:
        raise ValueError(f"--per-cluster must be between 1 and {PER_CLUSTER_DEFAULT}")
    blueprints = _build_blueprints()[:per_cluster]
    specs: list[CaseSpec] = []
    dataset_idx = 1
    for family in TASK_FAMILIES:
        for tier in COMPLEXITY_TIERS:
            for regime in SIGNAL_REGIMES:
                cluster_id = f"{family}__{tier}__{regime}"
                for bp_idx, bp in enumerate(blueprints, start=1):
                    paired_requested = family == "paired_coupled"
                    requested_mode = bp.base_mode
                    if paired_requested:
                        requested_mode = f"paired_coupled_{bp.base_mode}"
                    specs.append(
                        CaseSpec(
                            dataset_idx=dataset_idx,
                            dataset_id=f"gf_case_{dataset_idx:04d}",
                            cluster_id=cluster_id,
                            task_family=family,
                            complexity_tier=tier,
                            signal_regime=regime,
                            blueprint_index=bp_idx,
                            base_mode=bp.base_mode,
                            missing_level=bp.missing_level,
                            batch_level=bp.batch_level,
                            motif_driven=bp.motif_driven,
                            paired_requested=paired_requested,
                            requested_mode=requested_mode,
                            seed=BASE_SEED + dataset_idx,
                        )
                    )
                    dataset_idx += 1
    return specs


def _choose_from_range(rng: random.Random, values: list[int]) -> int:
    return values[rng.randrange(len(values))]


def _jitter(rng: random.Random, low: float, high: float) -> float:
    return round(rng.uniform(low, high), 6)


def _base_missing_range(tier: str, level: str) -> tuple[float, float]:
    table = {
        "simple": {"low": (0.0, 0.02), "medium": (0.08, 0.14), "high": (0.16, 0.22)},
        "moderate": {"low": (0.01, 0.04), "medium": (0.12, 0.20), "high": (0.22, 0.32)},
        "hard": {"low": (0.02, 0.06), "medium": (0.18, 0.28), "high": (0.30, 0.42)},
    }
    return table[tier][level]


def _base_batch_range(family: str, tier: str, level: str) -> tuple[float, float]:
    if family == "single_workflow":
        table = {
            "simple": {"light": (0.02, 0.08), "elevated": (0.45, 0.8)},
            "moderate": {"light": (0.03, 0.12), "elevated": (0.7, 1.2)},
            "hard": {"light": (0.05, 0.15), "elevated": (1.1, 2.0)},
        }
    elif family == "batch_correction_centric":
        table = {
            "simple": {"light": (0.45, 0.8), "elevated": (0.8, 1.4)},
            "moderate": {"light": (0.7, 1.2), "elevated": (1.1, 1.8)},
            "hard": {"light": (1.0, 1.6), "elevated": (1.6, 2.4)},
        }
    else:
        table = {
            "simple": {"light": (0.25, 0.55), "elevated": (0.6, 1.0)},
            "moderate": {"light": (0.45, 0.85), "elevated": (0.9, 1.4)},
            "hard": {"light": (0.7, 1.2), "elevated": (1.2, 2.0)},
        }
    return table[tier][level]


def _base_bio_range(regime: str) -> tuple[float, float]:
    table = {
        "absent_or_weak": (0.15, 0.75),
        "moderate": (0.9, 1.75),
        "strong_or_confounded": (1.8, 2.8),
    }
    return table[regime]


def _derive_generation_params(spec: CaseSpec) -> dict[str, Any]:
    rng = random.Random(spec.seed)
    glycan_choices = {
        "simple": [50, 80],
        "moderate": [80, 100],
        "hard": [100, 150],
    }
    sample_choices = {
        "simple": [20, 25, 30],
        "moderate": [15, 20, 25],
        "hard": [12, 15, 20],
    }
    n_glycans = _choose_from_range(rng, glycan_choices[spec.complexity_tier])
    n_samples_per_group = _choose_from_range(rng, sample_choices[spec.complexity_tier])
    missing_low, missing_high = _base_missing_range(spec.complexity_tier, spec.missing_level)
    batch_low, batch_high = _base_batch_range(spec.task_family, spec.complexity_tier, spec.batch_level)
    bio_low, bio_high = _base_bio_range(spec.signal_regime)

    missing_fraction = _jitter(rng, missing_low, missing_high)
    needs_missing = spec.missing_level != "low"
    if spec.task_family == "batch_correction_centric" and spec.complexity_tier == "simple":
        missing_fraction = _jitter(rng, 0.0, 0.04) if spec.missing_level == "low" else missing_fraction

    mnar_bias = _jitter(
        rng,
        1.05 if spec.missing_level == "low" else 1.35,
        1.4 if spec.missing_level == "low" else (2.6 if spec.complexity_tier != "hard" else 3.3),
    )

    kappa_mu = _jitter(rng, batch_low, batch_high)
    needs_batch = spec.task_family != "single_workflow" or spec.batch_level == "elevated"
    if spec.task_family == "single_workflow" and spec.batch_level == "light":
        needs_batch = False
        kappa_mu = _jitter(rng, 0.01, 0.08)

    var_b = _jitter(rng, max(0.05, batch_low * 0.55), max(0.12, batch_high * 0.95))
    bio_strength = _jitter(rng, bio_low, bio_high)
    confounded_signal = spec.signal_regime == "strong_or_confounded" and (
        spec.batch_level == "elevated" or spec.task_family != "single_workflow"
    )
    if confounded_signal:
        bio_strength = max(bio_strength, 1.9)
        kappa_mu = max(kappa_mu, 0.8)

    if spec.task_family == "batch_correction_centric":
        n_batches = 3 if spec.complexity_tier == "simple" else 4
    elif spec.task_family == "paired_coupled":
        n_batches = 2 if spec.complexity_tier == "simple" else 3
    else:
        n_batches = 2

    motif_mode = "motif_proxy" if spec.motif_driven else "feature_level"
    batch_mode = "confounded" if confounded_signal else ("required" if needs_batch else "light")
    missing_mode = "mnar_left_censored" if needs_missing else "minimal"
    return {
        "n_glycans": n_glycans,
        "n_H": n_samples_per_group,
        "n_U": n_samples_per_group,
        "n_batches": n_batches,
        "missing_fraction": missing_fraction,
        "mnar_bias": mnar_bias,
        "kappa_mu": kappa_mu,
        "var_b": var_b,
        "bio_strength": bio_strength,
        "needs_imputation": needs_missing,
        "needs_batch_correction": needs_batch,
        "confounded_signal": confounded_signal,
        "motif_mode": motif_mode,
        "batch_mode": batch_mode,
        "missing_mode": missing_mode,
    }


def _scenario_type(spec: CaseSpec, params: dict[str, Any]) -> str:
    if spec.task_family == "single_workflow":
        if params["needs_batch_correction"] and params["needs_imputation"]:
            return "single_workflow_mixed_mnar"
        if params["needs_batch_correction"]:
            return "single_workflow_batch_aware"
        if params["needs_imputation"]:
            return "single_workflow_missing_mnar"
        return "single_workflow_compositional_baseline"
    if spec.task_family == "batch_correction_centric":
        if params["confounded_signal"]:
            return "batch_correction_confounded"
        if params["needs_imputation"]:
            return "batch_correction_with_missing"
        return "batch_correction_focus"
    if params["confounded_signal"]:
        return "paired_coupled_confounded"
    if params["needs_imputation"] and params["needs_batch_correction"]:
        return "paired_coupled_mixed"
    return "paired_coupled_focus"


def _safe_unlink(path: Path) -> None:
    if path.exists():
        path.unlink()


def _extract_true_de_from_mask(meta: dict[str, Any], top_k: int = 20) -> list[str]:
    dp = meta.get("dirichlet_parameters", {}) or {}
    mask = dp.get("differential_mask") or dp.get("differential_masks")
    if isinstance(mask, list):
        ids = [f"G{i + 1:03d}" for i, v in enumerate(mask) if v]
        return ids[:top_k]
    return []


def _extract_true_de_from_alpha(meta: dict[str, Any], top_k: int = 20) -> list[str]:
    dp = meta.get("dirichlet_parameters", {}) or {}
    alpha_h = dp.get("alpha_H", []) or []
    alpha_u = dp.get("alpha_U", []) or []
    if not alpha_h or not alpha_u or len(alpha_h) != len(alpha_u):
        return []
    ranked = []
    for idx, (h, u) in enumerate(zip(alpha_h, alpha_u), start=1):
        ranked.append((abs(float(u) - float(h)), f"G{idx:03d}"))
    ranked.sort(reverse=True)
    return [feature_id for _, feature_id in ranked[:top_k]]


def _build_sample_metadata(meta: dict[str, Any], abundance_df: pd.DataFrame) -> pd.DataFrame:
    s_info = meta.get("sample_info", {}) or {}
    labels = s_info.get("bio_labels", []) or []
    batches = s_info.get("batch_labels", []) or []
    if "glycan_index" in abundance_df.columns:
        n = abundance_df.shape[1] - 1
    else:
        n = int(meta.get("n_samples", len(abundance_df)))
    if len(labels) != n:
        labels = [0] * (n // 2) + [1] * (n - n // 2)
    if len(batches) != n:
        n_batches = max(meta.get("batch_parameters", {}).get("n_batches", 2), 1)
        batches = [i % n_batches for i in range(n)]
    return pd.DataFrame(
        {
            "sample_id": [f"S{i + 1:03d}" for i in range(n)],
            "group_label": labels,
            "batch_label": batches,
        }
    )


def _mirror_optional_outputs(ds_dir: Path, seed: int) -> None:
    clean = ds_dir / f"1_Y_clean_seed{seed}.csv"
    if clean.exists():
        shutil.copy2(clean, ds_dir / "clean_abundance.csv")
    batch_a = ds_dir / f"2_Y_with_batch_seed{seed}.csv"
    if batch_a.exists():
        shutil.copy2(batch_a, ds_dir / "batch_abundance.csv")
    batch_c = ds_dir / f"2_Y_with_batch_clr_seed{seed}.csv"
    if batch_c.exists():
        shutil.copy2(batch_c, ds_dir / "batch_clr.csv")
    clean_clr = ds_dir / f"1_Y_clean_clr_seed{seed}.csv"
    if clean_clr.exists():
        shutil.copy2(clean_clr, ds_dir / "clean_clr.csv")
    for src in ds_dir.glob("*missing*.csv"):
        if "clr" in src.name:
            target = ds_dir / "missing_clr.csv"
        else:
            target = ds_dir / "missing_abundance.csv"
        if src.resolve() == target.resolve():
            continue
        shutil.copy2(src, target)


def _quality_flags(meta: dict[str, Any]) -> dict[str, Any]:
    qc = meta.get("quality_checks") or meta.get("quality") or {}
    return dict(qc) if isinstance(qc, dict) else {}


def _n_samples_n_features(abundance_df: pd.DataFrame) -> tuple[int, int, float, float]:
    if "glycan_index" in abundance_df.columns:
        numeric = abundance_df.drop(columns=["glycan_index"]).apply(pd.to_numeric, errors="coerce")
        n_samples = int(numeric.shape[1])
        n_features = int(numeric.shape[0])
        totals = numeric.fillna(0).sum(axis=0)
    else:
        numeric = abundance_df.apply(pd.to_numeric, errors="coerce")
        n_samples = int(numeric.shape[0])
        n_features = int(numeric.shape[1])
        totals = numeric.fillna(0).sum(axis=1)
    miss_frac = float(numeric.isna().sum().sum() / max(numeric.size, 1))
    comp_like = float(totals.std() / max(abs(totals.mean()), 1e-12)) if len(totals) else 1.0
    return n_samples, n_features, miss_frac, comp_like


def _build_feature_metadata(n_features: int, truth_features: list[str], motif_driven: bool) -> pd.DataFrame:
    motif_features = set(truth_features[: max(1, min(8, len(truth_features)))]) if motif_driven else set()
    motif_label = "alpha2-3_sialylation_proxy" if motif_driven else ""
    rows = []
    for idx in range(1, n_features + 1):
        glycan_id = f"G{idx:03d}"
        rows.append(
            {
                "glycan_id": glycan_id,
                "glycan_name": f"Glycan_{idx:03d}",
                "glycan_class": "N-glycan",
                "data_type": "relative_abundance",
                "motif_tag": motif_label if glycan_id in motif_features else "",
                "effect_role": "truth_signal" if glycan_id in truth_features else "background",
                "notes": "Stratified benchmark feature annotation",
            }
        )
    return pd.DataFrame(rows)


def _build_task_card(spec: CaseSpec, params: dict[str, Any], scenario_type: str) -> dict[str, Any]:
    family_text = {
        "single_workflow": "Recommend one coherent glycomics workflow from raw abundance matrix to interpretable DE results.",
        "batch_correction_centric": "Recommend a workflow with special attention to batch-effect diagnosis and correction without erasing biology.",
        "paired_coupled": "Recommend a workflow that respects paired or coupled sampling structure when possible and states fallback handling if pairing is imperfect.",
    }[spec.task_family]
    return {
        "task_id": spec.dataset_id,
        "task_type": "workflow_recommendation",
        "task_family": spec.task_family,
        "complexity_tier": spec.complexity_tier,
        "signal_regime": spec.signal_regime,
        "task_description": family_text,
        "analysis_goal": (
            "Identify differential glycan signals between two biological conditions while handling compositionality, "
            "missingness, and batch structure as needed."
        ),
        "data_files": {
            "abundance": "abundance.csv",
            "sample_metadata": "sample_metadata.csv",
            "feature_metadata": "feature_metadata.csv",
        },
        "known_constraints": {
            "data_type": "relative_abundance_glycomics",
            "compositionally_closed": True,
            "missing_values_present": bool(params["needs_imputation"]),
            "missingness_pattern": params["missing_mode"],
            "batch_labels_available": True,
            "batch_structure_expected": bool(params["needs_batch_correction"]),
            "paired_or_coupled_requested": bool(spec.paired_requested),
            "motif_driven_signal_requested": bool(spec.motif_driven),
            "scenario_type": scenario_type,
        },
        "schema": {
            "sample_id_col": "sample_id",
            "group_col": "group_label",
            "batch_col": "batch_label",
        },
        "expected_output": {
            "type": "structured_workflow_plan",
            "must_include": ["ordered_tool_sequence", "parameter_choices", "rationale"],
        },
    }


def _build_manifest(
    ds_dir: Path,
    spec: CaseSpec,
    params: dict[str, Any],
    scenario_type: str,
    task_card: dict[str, Any],
    sample_meta: pd.DataFrame,
    abundance_df: pd.DataFrame,
    tool_specs: list[dict[str, Any]],
    actual_mode: str,
    quality_flags: dict[str, Any],
    degradation: dict[str, Any] | None,
) -> dict[str, Any]:
    n_samples, n_features, miss_frac, comp_like = _n_samples_n_features(abundance_df)
    return {
        "dataset_id": spec.dataset_id,
        "task_type": "workflow_recommendation",
        "task_family": spec.task_family,
        "complexity_tier": spec.complexity_tier,
        "signal_regime": spec.signal_regime,
        "scenario_type": scenario_type,
        "paths": {
            "abundance": "abundance.csv",
            "clr": "clr.csv",
            "sample_metadata": "sample_metadata.csv",
            "clean_abundance": "clean_abundance.csv" if (ds_dir / "clean_abundance.csv").exists() else None,
            "clean_clr": "clean_clr.csv" if (ds_dir / "clean_clr.csv").exists() else None,
            "batch_abundance": "batch_abundance.csv" if (ds_dir / "batch_abundance.csv").exists() else None,
            "batch_clr": "batch_clr.csv" if (ds_dir / "batch_clr.csv").exists() else None,
            "missing_abundance": "missing_abundance.csv" if (ds_dir / "missing_abundance.csv").exists() else None,
            "missing_clr": "missing_clr.csv" if (ds_dir / "missing_clr.csv").exists() else None,
            "glycoforge_metadata": "glycoforge_metadata.json" if (ds_dir / "glycoforge_metadata.json").exists() else None,
            "task": "task.json",
            "feature_metadata": "feature_metadata.csv",
        },
        "schema": task_card["schema"],
        "dataset_summary": {
            "n_samples": n_samples,
            "n_features": n_features,
            "has_missing_values": miss_frac > 1e-9,
            "missing_fraction_observed": miss_frac,
            "has_batch_labels": True,
            "n_batches": int(sample_meta["batch_label"].nunique()),
            "is_compositional": bool(comp_like < 0.08 and n_features > 1),
            "has_clr_matrix": True,
            "source": "glycoforge",
        },
        "task": {
            "objective": task_card["analysis_goal"],
            "expected_output_type": "structured_plan_and_tool_trace",
            "user_question": (
                f"Recommend a glycomics workflow for task_family={spec.task_family}, "
                f"complexity_tier={spec.complexity_tier}, signal_regime={spec.signal_regime}, "
                f"scenario_type={scenario_type}."
            ),
        },
        "available_tools": tool_specs,
        "constraints": {
            "hidden_truth_forbidden": True,
            "max_steps": 20,
            "must_output_json_plan": True,
        },
        "properties": {
            "contains_abundance": True,
            "contains_clr": True,
            "contains_metadata": True,
            "contains_feature_metadata": True,
        },
        "generation": {
            "requested_mode": spec.requested_mode,
            "actual_mode": actual_mode,
            "base_mode": spec.base_mode,
            "params": {
                "n_glycans": params["n_glycans"],
                "n_H": params["n_H"],
                "n_U": params["n_U"],
                "n_batches": params["n_batches"],
                "missing_fraction": params["missing_fraction"],
                "mnar_bias": params["mnar_bias"],
                "kappa_mu": params["kappa_mu"],
                "var_b": params["var_b"],
                "bio_strength": params["bio_strength"],
            },
            "seed": spec.seed,
            "blueprint_index": spec.blueprint_index,
            "quality_flags": quality_flags,
            "degradation": degradation or {},
        },
    }


def _build_hidden_truth(
    spec: CaseSpec,
    params: dict[str, Any],
    actual_mode: str,
    sample_meta: pd.DataFrame,
    truth_features: list[str],
    degradation: dict[str, Any] | None,
) -> dict[str, Any]:
    needs_batch = bool(params["needs_batch_correction"])
    needs_missing = bool(params["needs_imputation"])
    recommended_batch_family = "run_combat" if needs_batch else "none"
    expected = ["include_clr_before_de_for_compositional_data"]
    expected.append("impute_when_missing" if needs_missing else "optional_impute")
    expected.append("batch_correct_before_de" if needs_batch else "optional_batch")
    if spec.motif_driven:
        expected.append("consider_motif_level_follow_up")
    if spec.paired_requested:
        expected.append("acknowledge_pairing_or_coupling")
    return {
        "dataset_id": spec.dataset_id,
        "task_family": spec.task_family,
        "complexity_tier": spec.complexity_tier,
        "signal_regime": spec.signal_regime,
        "true_conditions": {
            "batch_injected": needs_batch,
            "missingness_injected": needs_missing,
            "bio_strength": float(params["bio_strength"]),
            "batch_strength": float(params["kappa_mu"]),
            "missing_fraction_target": float(params["missing_fraction"]),
            "mnar_bias": float(params["mnar_bias"]),
            "motif_driven": bool(spec.motif_driven),
            "paired_or_coupled_requested": bool(spec.paired_requested),
            "paired_or_coupled_realized": _paired_realized(actual_mode, degradation),
            "confounded_signal": bool(params["confounded_signal"]),
        },
        "true_labels": {
            "bio_labels": sample_meta["group_label"].tolist(),
            "batch_labels": sample_meta["batch_label"].tolist(),
        },
        "true_differential": {
            "feature_ids": truth_features,
            "ranked": truth_features,
            "effect_level": "motif_proxy" if spec.motif_driven else "glycan",
        },
        "true_workflow_tags": {
            "needs_imputation": needs_missing,
            "needs_clr": True,
            "needs_batch_correction": needs_batch,
            "recommended_batch_family": recommended_batch_family,
            "should_preserve_bio_covariate": True,
        },
        "scoring_targets": {
            "expected_plan_requirements": expected,
            "disallowed_plan_patterns": [r"run_differential_expression.*run_combat"],
        },
        "generation": {
            "requested_mode": spec.requested_mode,
            "actual_mode": actual_mode,
            "degradation": degradation or {},
        },
        "glycoforge_metadata": {"path": "glycoforge_metadata.json"},
    }


def _copy_contract_files(src_dir: Path, dst_dir: Path) -> None:
    names = [
        "abundance.csv",
        "clr.csv",
        "sample_metadata.csv",
        "clean_abundance.csv",
        "clean_clr.csv",
        "batch_abundance.csv",
        "batch_clr.csv",
        "missing_abundance.csv",
        "missing_clr.csv",
        "glycoforge_metadata.json",
    ]
    for name in names:
        src = src_dir / name
        if src.exists():
            shutil.copy2(src, dst_dir / name)
        else:
            _safe_unlink(dst_dir / name)


def _load_source_pool(source_root: Path) -> list[SourceCase]:
    if not source_root.exists():
        return []
    pool: list[SourceCase] = []
    for dataset_dir in sorted(source_root.glob("gf_case_*")):
        manifest_path = dataset_dir / "manifest.json"
        truth_path = dataset_dir / "hidden_truth.json"
        if not (manifest_path.exists() and truth_path.exists()):
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        hidden_truth = json.loads(truth_path.read_text(encoding="utf-8"))
        task_path = dataset_dir / "task.json"
        task_card = json.loads(task_path.read_text(encoding="utf-8")) if task_path.exists() else None
        true_conditions = hidden_truth.get("true_conditions", {})
        generation = manifest.get("generation", {})
        dataset_summary = manifest.get("dataset_summary", {})
        pool.append(
            SourceCase(
                dataset_id=manifest.get("dataset_id", dataset_dir.name),
                dataset_dir=dataset_dir,
                manifest=manifest,
                hidden_truth=hidden_truth,
                task_card=task_card,
                missing_fraction=float(
                    true_conditions.get("missing_fraction_target", dataset_summary.get("missing_fraction_observed", 0.0))
                ),
                batch_strength=float(true_conditions.get("batch_strength", generation.get("params", {}).get("kappa_mu", 0.0))),
                bio_strength=float(true_conditions.get("bio_strength", generation.get("params", {}).get("bio_strength", 0.0))),
                source_mode=str(generation.get("actual_mode", generation.get("requested_mode", "synthetic"))),
                source_scenario=str(manifest.get("scenario_type", "")),
                has_missing=bool(true_conditions.get("missingness_injected", dataset_summary.get("has_missing_values", False))),
                has_batch=bool(true_conditions.get("batch_injected", False)),
            )
        )
    return pool


def _score_source_case(source: SourceCase, spec: CaseSpec, params: dict[str, Any]) -> tuple[float, str]:
    desired_missing = bool(params["needs_imputation"])
    desired_batch = bool(params["needs_batch_correction"])
    scenario_key = "clean"
    if desired_missing and desired_batch:
        scenario_key = "batch_and_missing"
    elif desired_missing:
        scenario_key = "missing_only"
    elif desired_batch:
        scenario_key = "batch_only"

    penalty = 0.0
    if desired_missing != source.has_missing:
        penalty += 2.0
    if desired_batch != source.has_batch:
        penalty += 2.0
    if spec.base_mode == "templated" and "templated" not in source.source_mode:
        penalty += 0.75
    if spec.base_mode == "synthetic" and "synthetic" not in source.source_mode:
        penalty += 0.25
    penalty += abs(source.missing_fraction - float(params["missing_fraction"])) * 4.0
    penalty += abs(source.batch_strength - float(params["kappa_mu"])) * 0.8
    penalty += abs(source.bio_strength - float(params["bio_strength"])) * 0.4
    if scenario_key == "batch_only" and source.source_scenario != "batch_only":
        penalty += 0.5
    if scenario_key == "missing_only" and source.source_scenario != "missing_only":
        penalty += 0.5
    if scenario_key == "batch_and_missing" and source.source_scenario != "batch_and_missing":
        penalty += 0.5
    return (penalty, source.dataset_id)


def _select_source_case(pool: list[SourceCase], spec: CaseSpec, params: dict[str, Any]) -> SourceCase:
    ranked = sorted(pool, key=lambda source: _score_source_case(source, spec, params))
    if not ranked:
        raise RuntimeError("No source cases available for clone fallback")
    offset = (spec.blueprint_index - 1) % min(len(ranked), 12)
    return ranked[offset]


def _detect_paired_kwargs(simulate_fn) -> dict[str, Any]:
    param_names = _simulate_param_names(simulate_fn)
    if "paired" in param_names:
        return {"paired": True}
    if "paired_design" in param_names:
        return {"paired_design": True}
    if "coupled" in param_names:
        return {"coupled": True}
    if "matched_pairs" in param_names:
        return {"matched_pairs": True}
    return {}


def _template_files(fallback_source_pack: Path) -> list[Path]:
    files: list[Path] = []
    files.extend(sorted(Path(".").glob("glycoforge_500_datasets/1_Y_clean_seed*.csv")))
    if fallback_source_pack.exists():
        files.extend(sorted(fallback_source_pack.glob("gf_case_*/clean_abundance.csv")))
    deduped: list[Path] = []
    seen: set[str] = set()
    for path in files:
        key = str(path.resolve())
        if key not in seen:
            seen.add(key)
            deduped.append(path)
    return deduped


def _paired_realized(actual_mode: str, degradation: dict[str, Any] | None) -> bool:
    if degradation and degradation.get("paired_coupled"):
        return False
    return actual_mode.startswith("paired_coupled_")


def _generate_with_simulate(
    spec: CaseSpec,
    params: dict[str, Any],
    ds_dir: Path,
    simulate_fn,
    simulate_paired_fn,
    template_files: list[Path],
) -> tuple[dict[str, Any], str, dict[str, Any] | None]:
    sim_kwargs: dict[str, Any] = {
        "data_source": "real" if spec.base_mode == "templated" else "simulated",
        "random_seeds": [spec.seed],
        "output_dir": str(ds_dir),
        "save_csv": True,
        "verbose": False,
        "n_glycans": params["n_glycans"],
        "n_H": params["n_H"],
        "n_U": params["n_U"],
        "n_batches": params["n_batches"],
        "missing_fraction": params["missing_fraction"],
        "mnar_bias": params["mnar_bias"],
        "kappa_mu": params["kappa_mu"],
        "var_b": params["var_b"],
        "bio_strength": params["bio_strength"],
    }
    degradation: dict[str, Any] = {}
    actual_mode = spec.base_mode

    if spec.base_mode == "templated" and template_files:
        sim_kwargs["use_real_effect_sizes"] = True
        sim_kwargs["data_file"] = str(template_files[(spec.dataset_idx - 1) % len(template_files)])
        sim_kwargs["column_prefix"] = {"healthy": "healthy_", "unhealthy": "unhealthy_"}

    if spec.paired_requested:
        if simulate_paired_fn is not None:
            degradation["paired_coupled"] = "simulate_paired_available_but_not_yet_materialized_in_contract"
            actual_mode = f"paired_proxy_from_{spec.base_mode}"
        else:
            degradation["paired_coupled"] = "requested_but_installed_glycoforge_has_no_simulate_paired"
            actual_mode = f"single_proxy_from_{spec.requested_mode}"

    try:
        simulate_fn(**sim_kwargs)
    except Exception as first_exc:
        if "use_real_effect_sizes" in sim_kwargs:
            degradation["templated"] = f"templated_failed:{type(first_exc).__name__}"
            sim_kwargs.pop("use_real_effect_sizes", None)
            sim_kwargs.pop("data_file", None)
            actual_mode = "synthetic"
            if spec.paired_requested and _detect_paired_kwargs(simulate_fn):
                actual_mode = "synthetic"
            simulate_fn(**sim_kwargs)
        else:
            raise

    metadata_path = ds_dir / f"metadata_seed{spec.seed}.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing metadata output: {metadata_path}")
    meta = json.loads(metadata_path.read_text(encoding="utf-8"))
    return meta, actual_mode, degradation or None


def _materialize_case_from_meta(
    spec: CaseSpec,
    params: dict[str, Any],
    ds_dir: Path,
    meta: dict[str, Any],
    actual_mode: str,
    tool_specs: list[dict[str, Any]],
    degradation: dict[str, Any] | None,
    source_case_id: str | None = None,
) -> dict[str, Any]:
    if params["needs_imputation"]:
        abundance_seed = ds_dir / "missing_abundance.csv"
        clr_seed = ds_dir / "missing_clr.csv"
        if not (abundance_seed.exists() and clr_seed.exists()):
            abundance_seed = ds_dir / f"3_Y_with_batch_and_missing_seed{spec.seed}.csv"
            clr_seed = ds_dir / f"3_Y_with_batch_and_missing_clr_seed{spec.seed}.csv"
    else:
        abundance_seed = ds_dir / "batch_abundance.csv"
        clr_seed = ds_dir / "batch_clr.csv"
        if not (abundance_seed.exists() and clr_seed.exists()):
            abundance_seed = ds_dir / f"2_Y_with_batch_seed{spec.seed}.csv"
            clr_seed = ds_dir / f"2_Y_with_batch_clr_seed{spec.seed}.csv"
    if abundance_seed.exists() and clr_seed.exists():
        abundance_df = pd.read_csv(abundance_seed)
        clr_df = pd.read_csv(clr_seed)
        abundance_df.to_csv(ds_dir / "abundance.csv", index=False)
        clr_df.to_csv(ds_dir / "clr.csv", index=False)
        metadata_seed = ds_dir / f"metadata_seed{spec.seed}.json"
        if metadata_seed.exists():
            shutil.copy2(metadata_seed, ds_dir / "glycoforge_metadata.json")
        _mirror_optional_outputs(ds_dir, spec.seed)
    else:
        abundance_df = pd.read_csv(ds_dir / "abundance.csv")

    sample_meta = _build_sample_metadata(meta, abundance_df)
    sample_meta.to_csv(ds_dir / "sample_metadata.csv", index=False)

    truth_features = _extract_true_de_from_mask(meta) or _extract_true_de_from_alpha(meta)
    scenario_type = _scenario_type(spec, params)
    feature_meta = _build_feature_metadata(
        n_features=_n_samples_n_features(abundance_df)[1],
        truth_features=truth_features,
        motif_driven=spec.motif_driven,
    )
    feature_meta.to_csv(ds_dir / "feature_metadata.csv", index=False)

    task_card = _build_task_card(spec, params, scenario_type)
    (ds_dir / "task.json").write_text(json.dumps(task_card, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = _build_manifest(
        ds_dir=ds_dir,
        spec=spec,
        params=params,
        scenario_type=scenario_type,
        task_card=task_card,
        sample_meta=sample_meta,
        abundance_df=abundance_df,
        tool_specs=tool_specs,
        actual_mode=actual_mode,
        quality_flags=_quality_flags(meta),
        degradation=degradation,
    )
    if source_case_id:
        manifest["generation"]["source_dataset_id"] = source_case_id
    (ds_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    hidden_truth = _build_hidden_truth(
        spec=spec,
        params=params,
        actual_mode=actual_mode,
        sample_meta=sample_meta,
        truth_features=truth_features,
        degradation=degradation,
    )
    if source_case_id:
        hidden_truth["generation"]["source_dataset_id"] = source_case_id
    (ds_dir / "hidden_truth.json").write_text(json.dumps(hidden_truth, ensure_ascii=False, indent=2), encoding="utf-8")

    return {
        "dataset_id": spec.dataset_id,
        "task_family": spec.task_family,
        "complexity_tier": spec.complexity_tier,
        "signal_regime": spec.signal_regime,
        "scenario_type": scenario_type,
        "task_type": "workflow_recommendation",
        "requested_mode": spec.requested_mode,
        "actual_mode": actual_mode,
        "seed": spec.seed,
        "n_glycans": params["n_glycans"],
        "n_samples": int(sample_meta.shape[0]),
        "missing_fraction": params["missing_fraction"],
        "kappa_mu": params["kappa_mu"],
        "var_b": params["var_b"],
        "bio_strength": params["bio_strength"],
        "motif_driven": spec.motif_driven,
        "paired_requested": spec.paired_requested,
        "paired_realized": _paired_realized(actual_mode, degradation),
        "source_dataset_id": source_case_id or "",
    }


def _run_case(
    spec: CaseSpec,
    tool_specs: list[dict[str, Any]],
    out_root: Path,
    simulate_fn,
    simulate_paired_fn,
    template_files: list[Path],
    source_pool: list[SourceCase],
    require_glycoforge: bool = False,
) -> dict[str, Any]:
    ds_dir = out_root / spec.dataset_id
    if ds_dir.exists():
        shutil.rmtree(ds_dir)
    ds_dir.mkdir(parents=True, exist_ok=True)
    params = _derive_generation_params(spec)

    if simulate_fn is not None:
        try:
            meta, actual_mode, degradation = _generate_with_simulate(
                spec=spec,
                params=params,
                ds_dir=ds_dir,
                simulate_fn=simulate_fn,
                simulate_paired_fn=simulate_paired_fn,
                template_files=template_files,
            )
            return _materialize_case_from_meta(
                spec=spec,
                params=params,
                ds_dir=ds_dir,
                meta=meta,
                actual_mode=actual_mode,
                tool_specs=tool_specs,
                degradation=degradation,
            )
        except Exception:
            if require_glycoforge or not source_pool:
                raise

    if not source_pool:
        raise RuntimeError("glycoforge is unavailable and no fallback source pack was found")

    source_case = _select_source_case(source_pool, spec, params)
    _copy_contract_files(source_case.dataset_dir, ds_dir)
    source_meta_path = ds_dir / "glycoforge_metadata.json"
    if source_meta_path.exists():
        meta = json.loads(source_meta_path.read_text(encoding="utf-8"))
    else:
        meta = {}
    degradation = {
        "strategy": "clone_existing_pack",
        "reason": "glycoforge_unavailable_or_generation_failed",
        "source_dataset_id": source_case.dataset_id,
    }
    actual_mode = "stratified_clone_existing_pack"
    if spec.paired_requested:
        actual_mode = "clone_proxy_for_paired_request"
        degradation["paired_coupled"] = "requested_but_represented_via_clone_proxy"
    return _materialize_case_from_meta(
        spec=spec,
        params=params,
        ds_dir=ds_dir,
        meta=meta,
        actual_mode=actual_mode,
        tool_specs=tool_specs,
        degradation=degradation,
        source_case_id=source_case.dataset_id,
    )


def _write_index(out_root: Path, rows: list[dict[str, Any]]) -> Path:
    index_path = out_root / "scenario_index.csv"
    fieldnames = [
        "dataset_id",
        "task_family",
        "complexity_tier",
        "signal_regime",
        "scenario_type",
        "task_type",
        "requested_mode",
        "actual_mode",
        "seed",
        "n_glycans",
        "n_samples",
        "missing_fraction",
        "kappa_mu",
        "var_b",
        "bio_strength",
        "motif_driven",
        "paired_requested",
        "paired_realized",
        "source_dataset_id",
    ]
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return index_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a 648-task stratified GlycoForge benchmark.")
    parser.add_argument("--out", type=str, default=str(DEFAULT_OUT_ROOT), help="Output directory.")
    parser.add_argument(
        "--per-cluster",
        type=int,
        default=PER_CLUSTER_DEFAULT,
        help="Tasks per (family, tier, regime) cluster. Default 24 -> 648 total.",
    )
    parser.add_argument(
        "--fallback-source-pack",
        type=str,
        default=str(DEFAULT_FALLBACK_SOURCE_PACK),
        help="Existing GlycoAgentBench dataset pack used for clone fallback when glycoforge is unavailable.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Keep the output directory if it exists; otherwise it is removed first.",
    )
    parser.add_argument(
        "--require-glycoforge",
        action="store_true",
        help="Fail instead of using clone fallback when glycoforge is unavailable or generation fails.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_root = Path(args.out)
    if out_root.exists() and not args.keep_existing:
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    specs = _build_case_specs(args.per_cluster)
    simulate_fn = _try_import_simulate()
    simulate_paired_fn = _try_import_simulate_paired()
    fallback_source_pack = Path(args.fallback_source_pack)
    template_files = _template_files(fallback_source_pack)
    source_pool = _load_source_pool(fallback_source_pack)
    tool_specs = export_tool_specs_for_manifest()

    if args.require_glycoforge and simulate_fn is None:
        raise RuntimeError(
            "glycoforge is not importable in the current Python environment. "
            "Use the project .venv interpreter or install glycoforge first."
        )

    if simulate_fn is None:
        print("INFO: glycoforge is not importable; using clone fallback when possible.")
    elif not template_files:
        print("INFO: templated mode requested for half the cases but no template files were found; synthetic fallback may be used.")
    if source_pool:
        print(f"INFO: fallback source pool ready: {len(source_pool)} cases from {args.fallback_source_pack}")

    rows: list[dict[str, Any]] = []
    for idx, spec in enumerate(specs, start=1):
        row = _run_case(
            spec=spec,
            tool_specs=tool_specs,
            out_root=out_root,
            simulate_fn=simulate_fn,
            simulate_paired_fn=simulate_paired_fn,
            template_files=template_files,
            source_pool=source_pool,
            require_glycoforge=args.require_glycoforge,
        )
        rows.append(row)
        if idx % max(1, min(54, len(specs))) == 0 or idx == len(specs):
            print(f"{idx}/{len(specs)} cases materialized")

    index_path = _write_index(out_root, rows)
    print(f"Done. Generated {len(rows)} stratified cases in {out_root}. Index: {index_path}")


if __name__ == "__main__":
    main()
