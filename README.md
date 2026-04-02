# GlycoForgeDataset

GlycoForgeDataset 是一个面向糖组学分析流程评测的合成数据集仓库。数据由 GlycoForge 生成或基于模板/回退机制物化，带有可验证的隐藏真值，适合用于评估计算方法、自动化分析流程和 LLM/Agent 在糖组学任务中的表现。

本仓库当前发布的数据包为 `glycoforge_benchmark_648/`，共 **648** 个 case，覆盖：

- 3 个任务家族：`single_workflow`、`batch_correction_centric`、`paired_coupled`
- 3 个复杂度等级：`simple`、`moderate`、`hard`
- 3 个信号强度层级：`absent_or_weak`、`moderate`、`strong_or_confounded`
- 每个 `(task_family, complexity_tier, signal_regime)` 组合下 24 个分层 case

## 1. 数据集目标

这个数据集主要用于以下场景：

- 评估糖组学分析流程推荐能力
- 评估缺失值处理、CLR 变换、批次校正和差异分析策略
- 评估模型或 Agent 是否能够根据数据属性选择合适工具链
- 在已知 ground truth 的条件下进行方法比较和误差分析

与普通公开数据不同，这个仓库中的每个 case 都不仅提供输入矩阵，还提供：

- 任务描述 `task.json`
- 机器可读清单 `manifest.json`
- 样本与特征注释
- 隐藏真值 `hidden_truth.json`，用于离线评测而不是模型输入

## 2. 数据集概览

发布包路径：

```text
glycoforge_benchmark_648/
```

其中包含：

- `scenario_index.csv`：全量 case 索引表
- `gf_case_0001/` 到 `gf_case_0648/`：每个 case 一个目录

按索引统计，当前 648 个 case 的分布是：

- 每个 `task_family` 各 216 个 case
- 每个 `complexity_tier` 各 216 个 case
- 每个 `signal_regime` 各 216 个 case

`actual_mode` 的分布如下：

- `synthetic`：216
- `templated`：216
- `degraded_from_paired_coupled_synthetic_to_synthetic`：108
- `degraded_from_paired_coupled_templated_to_templated`：108

这意味着数据集不仅覆盖标准合成生成，还显式记录了部分配对任务在当前生成能力下的退化实现状态。

## 3. 任务家族说明

### `single_workflow`

面向“给定一个糖组学相对丰度数据集，应如何选择一条合理分析流程”的任务。重点考察：

- 是否识别组合型数据特性
- 是否判断是否需要 CLR
- 是否识别缺失值问题
- 是否识别批次效应是否显著
- 是否能串联出完整、连贯的分析步骤

### `batch_correction_centric`

面向“批次效应处理是核心问题”的任务。重点考察：

- 是否优先检查并量化批次效应
- 是否选择合适的批次校正方法
- 是否避免过度校正
- 是否在校正后再进行下游差异分析

### `paired_coupled`

面向“样本存在配对/耦合结构需求”的任务。重点考察：

- 是否意识到配对设计对分析策略的影响
- 是否在配对无法完全实现时给出合理 fallback
- 是否在工作流建议中保留对耦合结构的约束意识

## 4. 每个 case 的目录结构

单个 case 目录示例：

```text
glycoforge_benchmark_648/
└── gf_case_0001/
    ├── abundance.csv
    ├── clr.csv
    ├── clean_abundance.csv
    ├── clean_clr.csv
    ├── batch_abundance.csv
    ├── batch_clr.csv
    ├── missing_abundance.csv
    ├── missing_clr.csv
    ├── sample_metadata.csv
    ├── feature_metadata.csv
    ├── glycoforge_metadata.json
    ├── task.json
    ├── manifest.json
    ├── hidden_truth.json
    ├── metadata_seed<seed>.json
    ├── 1_Y_clean_seed<seed>.csv
    ├── 1_Y_clean_clr_seed<seed>.csv
    ├── 2_Y_with_batch_seed<seed>.csv
    ├── 2_Y_with_batch_clr_seed<seed>.csv
    ├── 3_Y_with_batch_and_missing_seed<seed>.csv
    └── 3_Y_with_batch_and_missing_clr_seed<seed>.csv
```

主要文件含义如下：

| 文件 | 说明 |
| --- | --- |
| `abundance.csv` | 推荐作为主要输入的相对丰度矩阵 |
| `clr.csv` | 对应的 CLR 变换矩阵 |
| `sample_metadata.csv` | 样本注释，至少包含 `sample_id`、`group_label`、`batch_label` |
| `feature_metadata.csv` | 特征注释，包含 `glycan_id`、`glycan_name`、`motif_tag`、`effect_role` 等 |
| `task.json` | 当前 case 的任务卡片，描述目标、约束、预期输出 |
| `manifest.json` | 机器可读数据清单、统计摘要、工具列表、生成参数 |
| `glycoforge_metadata.json` | 生成过程元数据与质量检查结果 |
| `hidden_truth.json` | 真值标签与真实条件，适用于评测，不应暴露给模型 |
| `clean_*` / `batch_*` / `missing_*` | 对应不同退化阶段的数据快照 |
| `*_seed<seed>.csv` | 带随机种子的中间物化结果，便于复现与核查 |

## 5. 文件格式说明

### 5.1 `scenario_index.csv`

这是整个 benchmark 的总索引表，常见字段包括：

- `dataset_id`
- `task_family`
- `complexity_tier`
- `signal_regime`
- `scenario_type`
- `task_type`
- `requested_mode`
- `actual_mode`
- `seed`
- `n_glycans`
- `n_samples`
- `missing_fraction`
- `kappa_mu`
- `var_b`
- `bio_strength`
- `motif_driven`
- `paired_requested`
- `paired_realized`
- `source_dataset_id`

如果你要做批量实验，建议首先读取这个文件做筛选和分层抽样。

### 5.2 `sample_metadata.csv`

当前 case 中通常至少包含：

- `sample_id`：样本 ID
- `group_label`：生物学分组标签
- `batch_label`：批次标签

### 5.3 `feature_metadata.csv`

通常包含：

- `glycan_id`
- `glycan_name`
- `glycan_class`
- `data_type`
- `motif_tag`
- `effect_role`
- `notes`

其中 `effect_role` 可用于离线核对哪些特征属于真实信号。

### 5.4 `task.json`

用于给 Agent 或评测程序提供任务上下文，例如：

- 当前 case 的任务类型
- 数据文件入口
- 需要满足的输出结构
- 是否存在缺失值、批次、配对需求等先验约束

### 5.5 `manifest.json`

这是最完整的“对外说明书”，通常包含：

- 数据路径映射
- schema 说明
- 数据集摘要
- 任务目标
- 可用工具列表
- 生成参数
- 质量检查结果

如果你希望构建自动评测程序，`manifest.json` 往往是最合适的入口。

### 5.6 `hidden_truth.json`

该文件包含：

- 真实条件是否注入批次/缺失值
- 生物信号、批次信号、MNAR 参数
- 真实标签
- 真实差异特征或 motif 级信号

`hidden_truth.json` 只应在离线评估、打分或方法验证时使用，不应在模型推理阶段作为输入。

## 6. 推荐使用方式

### 方式 A：直接使用已发布数据

如果你的目标是：

- 跑 benchmark
- 评测工作流推荐
- 比较批次校正和差异分析策略
- 做 Agent 评测

那么直接使用仓库中的 `glycoforge_benchmark_648/` 即可，不需要重新生成。

最常见的入口文件是：

```text
glycoforge_benchmark_648/scenario_index.csv
```

以及某个具体 case 下的：

```text
abundance.csv
sample_metadata.csv
feature_metadata.csv
task.json
manifest.json
```

### 方式 B：重新生成 benchmark

仓库提供了生成脚本：

```text
generate_glycoforge_benchmark_648_stratified.py
```

但需要注意，这个脚本**不是完全自包含**的。它依赖：

- Python 环境
- `pandas`
- `glycoforge` 包
- 上级项目中的 `src.tools.registry`

也就是说，这个仓库更像是“数据发布仓库”，而不是一个独立可直接从零生成数据的最小工程。

## 7. 脚本运行说明

### 7.1 基本命令

如果你已经在完整的 GlycoForge 项目环境中，可以使用：

```bash
python generate_glycoforge_benchmark_648_stratified.py --out ./glycoforge_benchmark_648
```

说明：

- `--out` 指定输出目录
- 当前脚本默认输出目录是 `data/glycoforge_benchmark_648`
- 本仓库已发布的数据位于根目录下的 `./glycoforge_benchmark_648`
- 因此如果你希望复现当前仓库布局，建议显式传入 `--out ./glycoforge_benchmark_648`

### 7.2 常用参数

```bash
python generate_glycoforge_benchmark_648_stratified.py \
  --out ./glycoforge_benchmark_648 \
  --per-cluster 24
```

支持的参数包括：

- `--out`：输出目录
- `--per-cluster`：每个 `(task_family, complexity_tier, signal_regime)` 组合生成的 case 数量，默认 `24`，对应总数 `648`
- `--fallback-source-pack`：当 `glycoforge` 不可用或生成失败时，使用已有数据包作为回退来源
- `--keep-existing`：保留已有输出目录，不在运行前删除
- `--require-glycoforge`：要求必须成功导入并使用 `glycoforge`，否则直接报错，不允许回退

### 7.3 强制要求真实生成

如果你不希望脚本使用 clone/fallback 模式，而是要求必须调用可用的 `glycoforge`：

```bash
python generate_glycoforge_benchmark_648_stratified.py \
  --out ./glycoforge_benchmark_648 \
  --require-glycoforge
```

### 7.4 使用 fallback source pack

当 `glycoforge` 不可用时，脚本会尝试从已有数据包中复制并重建 contract 文件。示例：

```bash
python generate_glycoforge_benchmark_648_stratified.py \
  --out ./glycoforge_benchmark_648 \
  --fallback-source-pack /path/to/data/glycoforge_benchmark_500
```

## 8. 生成逻辑摘要

脚本使用分层设计组合生成 `CaseSpec`，主要维度包括：

- `task_family`
- `complexity_tier`
- `signal_regime`
- `base_mode`
- `missing_level`
- `batch_level`
- `motif_driven`

在生成过程中，会为每个 case 派生如下核心参数：

- `n_glycans`
- `n_samples_per_group`
- `missing_fraction`
- `mnar_bias`
- `kappa_mu`
- `var_b`
- `bio_strength`
- `n_batches`

随后脚本会物化：

- 主输入矩阵
- CLR 矩阵
- clean/batch/missing 阶段矩阵
- 样本和特征注释
- 任务卡片
- manifest
- hidden truth
- 总索引文件

## 9. 分析流程建议

如果你是首次使用这个数据集，推荐按如下顺序：

1. 读取 `scenario_index.csv`，选择一个子集
2. 对每个 case 读取 `task.json` 和 `manifest.json`
3. 以 `abundance.csv` 和 `sample_metadata.csv` 作为主要分析输入
4. 视情况使用 `clr.csv` 或自行从 abundance 进行 CLR
5. 完成 workflow recommendation 或下游分析
6. 最后再使用 `hidden_truth.json` 做离线打分

一个最简单的实验单元通常只需要：

- `abundance.csv`
- `sample_metadata.csv`
- `feature_metadata.csv`
- `task.json`

## 10. 使用注意事项

- 这是合成数据集，适合做方法评测与流程验证，不应直接替代真实生物学结论
- `hidden_truth.json` 不应暴露给待评测模型
- `paired_coupled` 任务中，`paired_requested` 与 `paired_realized` 可能不同，建议结合 `scenario_index.csv` 和 `manifest.json` 解读
- 某些 case 的 `actual_mode` 可能是退化实现，使用时应记录这一点
- 生成脚本依赖外部项目环境；如果只想使用数据，不需要运行脚本

## 11. 仓库结构

```text
.
├── glycoforge_benchmark_648/
│   ├── scenario_index.csv
│   ├── gf_case_0001/
│   ├── gf_case_0002/
│   └── ...
├── generate_glycoforge_benchmark_648_stratified.py
├── LICENSE
└── README.md
```

## 12. License

本仓库使用 [MIT License](./LICENSE)。

## 13. 引用建议

如果你在论文、报告或 benchmark 中使用本数据集，建议至少说明以下信息：

- 数据集名称：`GlycoForgeDataset`
- 发布子集：`glycoforge_benchmark_648`
- 使用的 case 筛选条件
- 是否使用 `hidden_truth.json` 进行离线评估
- 是否排除了 `paired_requested != paired_realized` 的 case

如果后续你希望，我还可以继续补两类内容：

1. 一个更偏 GitHub 展示风格的英文版 README
2. 一个更偏论文/benchmark 附录风格的数据字典章节
