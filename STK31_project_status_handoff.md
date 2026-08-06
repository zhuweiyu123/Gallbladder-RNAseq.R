# STK31 胆囊癌单细胞项目 — 现状交接包

> 用途：把当前项目进度、已完成内容、关键缺口与下一步建议打包，便于交给 GPT / 其他协作者继续讨论或修改。
>
> 生成日期：2026-07-17  
> 修订：2026-07-17（Codex 复核 + refined T/NK 冻结 + 样本级验证收尾）
>
> 项目路径（本地）：`E:\ZHUWEIYU\Documents\R`
>
> 运行环境（服务器）：`/home/zhuweiyu/codex-r`（见 `README_codex.md`）

---

## 0. 当前统一结论（冻结）

> 在严格拆分 T/NK 后，高置信 NK 从 7,419 降至 1,551。refined 分析仍提示 STK31-high 肿瘤上皮与 NK 之间存在候选互作轴，但当前 MHC-I、TIGIT/NECTIN、TGF-beta 和 NKG2D 的 high-vs-low 差异均为 similar，尚无证据支持 STK31-high 特异增强，也不支持因果关系。

### 0.1 Refined 标签已冻结（P0 已完成）

| 项 | 值 |
|----|-----|
| 正式字段 | **`analysis_celltype`** |
| 正式 NK | `analysis_celltype == "NK_cell"` |
| refined 高置信 NK | **1,551** |
| 旧 NK（cluster 0/1） | **7,419** |
| 保留比例 | **~20.8%** |
| 交集 / 旧有新无 / 新有旧无 | 1,546 / 5,873 / 5 |
| 冻结 RDS | `results/merged_tnk_refined_annotation/refined_merged_object.rds` |
| 冻结说明 | `results/merged_tnk_refined_annotation/refined_label_freeze_note.md` |

**cluster 0/1 并非纯 NK**（旧注释明显过宽）：

| 原 cluster | NK | T | ambiguous | contaminant |
|------------|-----|-----|-----------|-------------|
| 0 | 877 | 1,602 | 994 | 375 |
| 1 | 669 | 2,274 | 329 | 299 |

**禁止**：用 cluster 0/1 或旧 `manual_celltype=="NK_cell"` 定义正式 NK。

### 0.2 Refined 主结果目录（current）

- `results/merged_tnk_refined_annotation/`
- `results/merged_stk31_nk_refined_analysis/`
- `results/merged_cellchat_refined_high_low_nk/`
- `results/refined_sample_level_validation_v2/`（样本级验证 current；探索性）
- 清单：`results/STK31_refined_result_manifest.csv`（legacy / superseded 不得用于主结论）
- `results/refined_sample_level_validation/` → **superseded_method_issue**（v1 方法问题，仅审计）

基于**旧 NK** 的功能评分、DE、CellChat → 状态 **legacy**，以降级；以 refined 为准。

### 0.3 STK31 与机制轴（refined CellChat）

- STK31-high 上皮仍为 **116**，阳性率极低（三种定义敏感性一致）。
- MHC-I：high→NK=**0.1424**，low→NK=**0.1331**，delta=**+0.0093**，`similar`
- TIGIT/NECTIN、TGFb、NKG2D：均为 **`similar`**
- IFNG：`not_detected`
- **不能**说 STK31-high 特异增强；**不能**说因果调控。

### 0.4 已完成、勿再当缺口

- 上皮 STK31 **high vs low DE/GO**（已完成；不依赖 NK 标签）
- CellChat 含 high/low/NK 及 **正式 high-vs-low 差分表/图**（refined 目录）
- 观察性方向 bootstrap / LOSO（虚拟 KO 与样本级验证中均有；仍非因果）
- T/NK 严格拆分与标签冻结

### 0.5 样本级验证 v2（script 07 方法修订后，探索性）

目录：`results/refined_sample_level_validation_v2/`（单位=sample，n=4）  
v1 `refined_sample_level_validation/` 已标 **superseded_method_issue**（方向错误、检测定义、bootstrap/LOSO 口径问题）。

| 轴 | 方向 | 四样本标签 | direction_consistent | loso_majority_pattern_stable | within_sample_bootstrap_consistency | evidence |
|----|------|------------|----------------------|------------------------------|-------------------------------------|----------|
| MHC-I/HLA_epi_to_NK | epi→NK | 3 similar, 1 higher_in_low | FALSE | TRUE | some_samples_stable | exploratory_mixed |
| NECTIN2/PVR_epi_to_NK | epi→NK | 2 high, 1 similar, 1 low | FALSE | FALSE | some_samples_stable | exploratory_mixed |
| MICA/B_epi_to_NK | epi→NK | 2 high, 1 similar, 1 low | FALSE | FALSE | all_evaluable_samples_stable | exploratory_mixed |
| IFNG_NK_to_epi | NK→epi | 3 high, 1 low | FALSE | TRUE | some_samples_stable | exploratory_mixed |
| TGFb_NK_to_epi | NK→epi | 2 high, 2 similar | FALSE | FALSE | some_samples_stable | exploratory_mixed |
| TGFb_epi_to_NK | epi→NK | 2 high, 2 similar | FALSE | FALSE | some_samples_stable | exploratory_mixed |

要点：
- **不得**用 LOSO stable 覆盖 `direction_consistent=FALSE`
- within-sample bootstrap **不得**升级 evidence_level；非样本级推断
- STK31 三种定义细胞集合完全一致 → 敏感性 **non-informative**
- 主结论仍以 refined CellChat 为准：主轴 similar 候选；IFNG CellChat not_detected；无特异增强、无因果
- 样本级 IFNG_NK_to_epi 为 mixed（3 high/1 low），是 IFNGR 在 high/low 上皮差异的探索性代理，**不升格**主结论

### 0.6 仍存短板

- 批次未整合；恶性上皮未确认；外部验证不足；n=4 样本级仅为探索性
- scTenifoldKnk 网络 outdegree=0，不可作因果 KO

---

## 0b. 项目一句话

胆囊癌 10x scRNA-seq（tissue1/2/4/5）：**STK31-high 上皮 ↔ refined NK** 的**候选**互作轴（HLA-I / TIGIT-NECTIN / TGF-β / NKG2D）。  
T/NK 已严格拆分并冻结；机制证据仍为候选级（high-vs-low 多为 similar）。

---

## 1. 数据与仓库结构

### 1.1 数据

| 项目 | 内容 |
|------|------|
| 疾病 | 胆囊癌 gallbladder cancer |
| 样本 | tissue1, tissue2, tissue4, tissue5 |
| 格式 | 10x `filtered_feature_bc_matrix` |
| 本地数据目录 | `胆囊癌/`（或服务器 `gallbladder_cancer` 软链） |
| 服务器工作目录 | `/home/zhuweiyu/codex-r` |

### 1.2 脚本（按执行顺序）

| 脚本 | 作用 | 主要输出 |
|------|------|----------|
| `scripts/01_tissue2_basic_seurat.R` | tissue2 单样本基础 Seurat | `results/tissue2_basic_seurat/` |
| `scripts/02_merged_basic_analysis.R` | 四样本合并 + STK31/NK 主分析 + CellChat/GO | `merged_basic_seurat/`、`merged_stk31_nk_analysis/`、CellChat 相关目录 |
| `scripts/03_tissue2_stk31_nk_analysis.R` | tissue2 STK31/NK 独立复现 + validation | `tissue2_stk31_nk_analysis/` |
| `scripts/04_create_stk31_split_umap.R` | STK31 split UMAP 补图 | 相关 UMAP |
| `scripts/05_stk31_virtual_knockout.R` | scTenifoldKnk 虚拟敲除（观察性；网络不可解释） | `merged_stk31_virtual_knockout/` |
| `scripts/06_refine_t_nk_and_rerun_stk31_nk.R` | **T/NK 严格拆分 + refined DE/CellChat** | `merged_tnk_refined_annotation/` 等 |
| `scripts/07_refined_sample_level_validation.R` | 样本级验证 v2（方向感知/检测阈值） | `refined_sample_level_validation_v2/` |
| `scripts/create_stk31_nk_mentor_ppt.py` | 旧导师 PPT（legacy） | `STK31_NK_mentor_ppt/` |
| `scripts/create_stk31_nk_mentor_ppt_refined.py` | **refined 导师 PPT** | `STK31_NK_导师汇报_refined版.pptx` |

### 1.3 文档

- `README.md`：脚本顺序与数据策略（不提交 raw / 大 rds）
- `README_codex.md`：SSH 服务器、路径、运行规范
- `AGENTS.md`：token 优化与修改规范（工程助手约束）

---

## 2. 当前分析流程（端到端）

```text
原始 10x 矩阵 (tissue1/2/4/5)
        │
        ├─[01] tissue2 基础 Seurat
        │       QC → Normalize → HVGs → PCA → Cluster → UMAP
        │
        ├─[02] 四样本 merge 主分析（当前主战场）
        │       ├─ 合并 QC / 降维聚类（直接 merge，无 Harmony/CCA）
        │       ├─ STK31 定位（上皮内 high/low）
        │       ├─ NK 候选 cluster 鉴定（marker score）
        │       ├─ 粗注释 + 手动注释
        │       ├─ DE（STK31 high vs low 上皮；high 上皮 vs NK/all；NK vs all）
        │       ├─ GO BP 富集
        │       ├─ 候选 ligand-receptor 初筛（自定义参考表）
        │       ├─ follow-up：样本级关系、NK 功能模块评分
        │       └─ CellChat：全量 + STK31-high epithelial ↔ NK 聚焦
        │
        ├─[03] tissue2 平行机制线 + validation/
        │
        ├─[04] 可视化补图
        │
        ├─[05] scTenifoldKnk 虚拟敲除（上皮内 STK31）
        │
        └─ 导师汇报 PPT 打包
```

### 2.1 关键参数（基础流程）

- min.features 200；max.features 6000；max percent.mt 20
- 2000 HVGs；30 PCs；用 1:20 dims；resolution 0.5
- STK31-high：上皮内分位数阈值（summary 中 cutoff 显示为 0，即按检出/阈值逻辑定义 high）
- NK：marker score 自动选 cluster（merged 上为 0 和 1）

---

## 3. 已完成内容（有结果文件支撑）

### 3.1 基础单细胞

- 合并后：**28,751 cells**，**19 clusters**，**32,535 genes**
- QC 表、UMAP（sample / cluster）、cluster × sample 计数
- tissue2 单样本线已独立跑通

### 3.2 STK31 定位与注释

| 指标 | Merged 结果 |
|------|-------------|
| 上皮细胞 | 4,236 |
| STK31-high / positive 上皮 | **116（约 2.7%）** |
| STK31-low 上皮 | 4,120 |
| NK 候选细胞 | **7,419**（cluster 0,1） |
| 其他细胞 | 17,096 |

样本级 STK31 阳性（全细胞 pct）：

| sample | pct_STK31+ | 上皮 STK31+ pct | NK 占比 |
|--------|------------|-----------------|--------|
| tissue1 | ~0.73% | ~3.7% | ~29% |
| tissue2 | ~1.19% | ~3.5% | ~18% |
| tissue4 | ~0.85% | ~0.85% | ~16% |
| tissue5 | ~0.68% | ~10.1%（上皮总数少） | ~44% |

注释：Epithelial / NK / T / Myeloid / Macrophage / B / Plasma / Fibroblast / Endothelial / Mast 等已有 CSV + UMAP。

### 3.3 DE / GO / 功能评分（含 high vs low — 已完成）

**已完成**（Grok 初评曾弱化此项，Codex 纠正）：

- `markers_stk31_high_vs_low_tumor_epithelial.csv` ← **上皮内 STK31 high vs low DE**
- `stk31_high_vs_low_tumor_epithelial_go_bp.csv` ← **对应 GO**
- `markers_stk31_high_tumor_epithelial_vs_nk.csv`
- `markers_stk31_high_vs_all_other.csv`
- `markers_nk_vs_all_other.csv`
- 其他 GO BP 表
- NK 功能模块：cytotoxicity / activation / checkpoint / migration / TGFbeta_response
- 样本级 STK31–NK 关系表

### 3.4 细胞互作

#### A. 自定义候选 L-R 得分（靠前，仅作候选）

| 方向 | 通路 | 表述边界 |
|------|------|----------|
| NK → STK31-high epi | TGF beta suppression | 候选轴 |
| STK31-high epi → NK | TIGIT checkpoint | 候选轴 |
| NK → epi | IFNG response | 候选轴 |
| STK31-high epi → NK | NKG2D stress ligand | 候选轴 |
| STK31-high epi → NK | MHC I inhibitory KIR | **非 STK31-high 特异**（见下） |

#### B. CellChat 实际状态（Codex 精确表述）

- CellChat **已同时纳入** STK31-high 上皮、STK31-low 上皮、NK 等分组并计算通讯
- **尚未正式输出**：high-vs-low **差异互作统计**与对比图
- 因此缺口不是“完全没算 low”，而是“**没做正式的 high vs low 差分互作报告**”

聚焦机制轴（候选）：

- MHC-I / HLA
- TIGIT / NECTIN
- TGFb
- IFNG / IFN

**MHC-I 特别警告（Codex）**：

> high→NK 的 MHC-I 总概率并不高于 low→NK，**暂时不能**表述为“STK31-high 特异增强 MHC-I 抑制”。

整体结论上限：

> **仅支持“存在候选互作轴”**，不支持特异增强或因果。

### 3.5 虚拟敲除 + 观察性方向分析

`results/merged_stk31_virtual_knockout/analysis_summary.csv`：

| 项 | 值 |
|----|-----|
| epithelial_cells | 4236 |
| stk31_positive_cells | 116 |
| network_genes | 1500 |
| **stk31_network_outdegree** | **0** |
| **stk31_connected_genes** | **0** |
| network_interpretation | **not_interpretable** |
| direction_interpretation | Association-based；需实验验证 |

分层理解：

1. **网络敲除本体**：不可解释（outdegree=0）→ 不能当因果 KO
2. **观察性方向分析（已做，Grok 初评漏写）**：
   - 按样本匹配
   - bootstrap
   - leave-one-sample-out
   - 相关输出：`stk31_ko_sample_stability.csv`、`stk31_ko_leave_one_sample_out.csv`、方向汇总图等
3. 机制基因多数 `uncertain`；priority 多标 `not_interpretable_zero_STK31_connectivity`

### 3.6 汇报材料

- `results/STK31_NK_mentor_ppt/STK31_NK_导师汇报_图选择与机制解释.pptx` 已生成

---

## 4. 缺少的关键分析（按重要性，Codex 修订版）

### 4.1 P0 — 必须先做

1. **严格拆分 T / NK（最高优先）**  
   - 证据：`cluster_celltype_annotation.csv` 中 cluster 0、1  
     - `broad_celltype` = **T_cell**  
     - `manual_celltype` / `cluster_celltype_label` = **NK_cell**（直接覆盖）  
   - 影响：NK 功能评分、DE、CellChat、汇报图全部下游  
   - 动作：CD3D/CD3E/TRAC vs KLRD1/NCAM1/FCGR3A/NKG7 等 DotPlot + 必要时亚聚类后再标注

### 4.2 P1 — 身份过关后立刻做

2. **CellChat high-vs-low 差异互作正式输出**  
   - 数据已在模型里（high / low / NK 都算过）  
   - 缺：差分统计表 + 对比图 + 明确“特异 vs 上皮通用”

3. **STK31 阳性定义敏感性分析**  
   - 阳性率极低；counts≥1 / 分位数 / 按样本分别定义

4. **表述收敛**  
   - 只写“候选互作轴”  
   - **禁止**在现有数字下写“STK31-high 特异增强 MHC-I 抑制”

### 4.3 P2 — 增强可信度

5. 批次整合（Harmony/RPCA）后重注释、重跑关键结论  
6. 恶性上皮确认（inferCNV 等）  
7. 外部验证（bulk/公共 scRNA/预后）  
8. 实验设计（IF / 共培养 / KO-OE）  
9. 样本级统计加强（在已有 bootstrap / LOSO 基础上，可加伪 bulk 等）

### 4.4 已完成、勿再当缺口

- 上皮 STK31 **high vs low DE/GO**（已有）
- CellChat 中 **同时计算 high 与 low**（已有计算，缺差分报告）
- 观察性方向的 **样本匹配 / bootstrap / leave-one-sample-out**（已有）

### 4.5 可选

- Doublet 去除、系统通路评分、空间组学、velocity 等

---

## 5. 下一步建议（可执行，修订版）

### 5.1 唯一正确的下一步顺序

```text
P0  严格 T/NK 拆分与重注释
    → 冻结一套 celltype labels

P1  用新 labels 重跑：
    - NK 功能模块
    - 与 NK 相关的 DE
    - CellChat（含 high / low / NK）
    - 正式输出 high-vs-low 差异互作表+图

P1  同步：STK31 定义敏感性 + 证据等级表
    - MHC-I：对比 high→NK vs low→NK，禁止“特异增强”除非数字支持

P2  批次整合 / inferCNV / 外部验证 / 实验
```

### 5.2 机制主轴（讨论用，均为候选）

1. MHC-I / HLA → KIR（**注意：当前未见 high > low**）
2. NECTIN2/PVR → TIGIT
3. TGF-β
4. MICA/B → NKG2D

---

## 6. 给后续模型的协作指令（可直接复制）

```text
基于 STK31_project_status_handoff.md（含 Codex 复核修订）继续工作。

已共识：
1. 不要重做基础 Seurat 流程。
2. scTenifoldKnk 网络不可解释；方向分析仅观察性（已有样本匹配/bootstrap/LOSO）。
3. 上皮 STK31 high vs low 的 DE/GO 已完成，不要再说“缺 high vs low DE”。
4. CellChat 已计算 high/low/NK；缺的是正式 high-vs-low 差异互作统计与对比图。
5. MHC-I：当前不能说 STK31-high 特异增强（high→NK 不高于 low→NK）。
6. 结论上限：仅“存在候选互作轴”。

当前唯一 P0：
- 严格拆分 T/NK。cluster 0/1 broad=T_cell，manual 被直接标成 NK_cell。
- 在新注释冻结前，不扩大 NK 功能/CellChat 主结论。

P0 之后：
- 重跑 NK 相关模块与 CellChat
- 输出 high-vs-low 差异互作表+图
- 再考虑批次整合 / inferCNV / 外部验证

改代码：优先新增 scripts/06_xxx.R；服务器路径 /home/zhuweiyu/codex-r；结果写 results/ 新目录，日志写 logs/。
```

---

## 7. 关键结果文件速查

```text
results/
├── merged_basic_seurat/
│   ├── analysis_summary.csv
│   └── gallbladder_cancer_merged_basic_seurat.rds
├── merged_stk31_nk_analysis/
│   ├── analysis_summary.csv
│   ├── cluster_celltype_annotation.csv
│   ├── candidate_pathway_summary.csv
│   ├── markers_*.csv
│   ├── *go_bp.csv
│   └── followup/
├── merged_stk31_nk_followup/
│   ├── followup_summary.csv
│   ├── sample_level_stk31_nk_relationship.csv
│   └── nk_function_scores_by_sample.csv
├── merged_cellchat_focused_stk31_nk/
│   └── stk31_high_epithelial_nk_mechanism_axis_summary.csv
├── merged_cellchat_go_plots/
│   └── merged_cellchat_object.rds
├── merged_stk31_virtual_knockout/
│   ├── analysis_summary.csv   # outdegree=0，关键
│   └── stk31_ko_priority_candidates.csv
├── tissue2_basic_seurat/
├── tissue2_stk31_nk_analysis/
├── tissue2_stk31_nk_validation/
└── STK31_NK_mentor_ppt/
```

---

## 8. 当前结论的证据等级（Codex 对齐版）

| 结论 | 证据等级 | 说明 |
|------|----------|------|
| STK31 在少数上皮细胞中检出 | 中 | 有表与图；阳性率极低 |
| 上皮内 STK31 high vs low 转录差异 | 中 | **DE/GO 已完成** |
| cluster 0/1 是“纯 NK” | **低** | broad=T_cell，被手动覆盖为 NK；**未过关** |
| 存在候选互作轴（TIGIT/MHC-I/TGF-β/NKG2D 等） | 中偏低 | L-R + CellChat；**仅候选** |
| STK31-high **特异**增强 MHC-I→NK 抑制 | **不支持** | high→NK 不高于 low→NK |
| STK31 因果调控上述通路 | **低** | 网络 KO 不可解释；观察性方向有 bootstrap/LOSO 仍非因果 |
| 临床靶点 | **不足** | 无外部与功能实验 |

---

## 9. 一句话总结（修订）

计算主线已齐（含 high vs low DE/GO、CellChat 含 low 组、观察性方向的样本稳健性），但 **cluster 0/1 的 T/NK 身份未严格拆分**，这是当前最大阻塞；在此之前，结果只支持 **“存在候选互作轴”**，且 **不能** 声称 STK31-high 特异增强 MHC-I 抑制。下一步唯一正确顺序：**先 T/NK 拆分 → 再重跑 NK/CellChat 并输出 high-vs-low 差异互作报告**。
```
