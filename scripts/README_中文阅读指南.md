# 脚本中文阅读指南

这个目录包含 21 个脚本：17 个 R 分析脚本、3 个 Python 绘图/汇报脚本和 1 个 PowerShell 图片工具。每个脚本开头都有用途、输入、流程和输出说明，关键函数与步骤旁补充了中文注释。

建议先看懂 **01 → 03**，再阅读多样本流程 **02 → 06 → 07**。这些脚本包含不同阶段的分析版本，**不要把 01–18 理解成必须从头到尾依次运行的唯一流程**。

## 每个脚本做什么

下面的输入只是主要依赖；实际运行还要满足脚本顶部路径、对象字段、所需 R 包和外部数据的要求。

| 脚本 | 主要用途 | 主要前置结果 |
|---|---|---|
| [01](01_tissue2_basic_seurat.R) | tissue2 质控、聚类和 UMAP | tissue2 原始 10x 数据 |
| [02](02_merged_basic_analysis.R) | 四样本合并、STK31/NK 分析、通讯及后续图 | 四个样本原始 10x 数据 |
| [03](03_tissue2_stk31_nk_analysis.R) | tissue2 的 STK31/NK 与身份验证 | 01 的基础对象 |
| [04](04_create_stk31_split_umap.R) | 使用已有 UMAP 画 STK31 高低上皮 | 02/03 的对象、注释与汇总 |
| [05](05_stk31_virtual_knockout.R) | 虚拟敲除及观察性方向探索 | 02 的对象和细胞注释 |
| [06](06_refine_t_nk_and_rerun_stk31_nk.R) | 细化 T/NK 标签并重跑下游 | 02 的对象和旧注释 |
| [07](07_refined_sample_level_validation.R) | 样本级方向、留一和重采样验证 | 06 的冻结注释对象 |
| [08](08_plot_refined_candidate_lr_bubble.R) | 人工候选配体受体气泡图 | 06 的候选配体受体表 |
| [09](09_refined_nk_go_enrichment.R) | NK 差异基因 GO | 06 的 NK 差异表 |
| [10](10_refined_stk31_high_vs_low_and_all_go.R) | 补充 high 上皮差异与 GO | 06 的冻结注释对象 |
| [11](11_integrate_cnv_binary_and_rerun_core.R) | CNV 注释整合、差异、患者配对伪批量 | 外部新对象/CNV 表及 06 的旧注释 |
| [12](12_validate_cnv_binary_integration.R) | 独立回读验收 | 11 保存的对象与结果 |
| [13](13_make_p1_p4_requested_figures.R) | P1/P4 上皮 STK31/HLA 图 | 11 的整合对象 |
| [14](14_update_p1_p2_p4_and_nk_receptor_figures.R) | P1/P2/P4 上皮与本地 NK 受体图 | 11 的整合对象和保留的旧 NK 标签 |
| [15](15_update_public_qc4_hlae_kir2dl2_figures.R) | 公共 QC4 的 HLA-E/KIR2DL2 图 | 公共项目目标矩阵、元数据及 NK 对象 |
| [16](16_qc4_and_sample_cellchat_pathways.R) | 公共 QC4/本地各样本通讯及 WNT | 公共 10x/提取程序及 11 的对象 |
| [17](17_make_wnt_stk31_hla_figures.py) | 从表格制作 WNT/HLA 专题图 | 16、07、11 的相应 CSV |
| [18](18_run_missing_cnv_binary_de_go_cellchat.R) | CNV 分组表达图、GO/GSEA、CellChat | 11 的对象和差异表 |
| [旧版 PPT](create_stk31_nk_mentor_ppt.py) | 将旧结果图表拼成汇报 | 本机历史图表与摘要 |
| [refined PPT](create_stk31_nk_mentor_ppt_refined.py) | 将细化注释结果拼成汇报 | 06/07 等产生的图表 |
| [图片工具](imagine.ps1) | 调用接口生成并下载图片 | 图片描述与接口配置 |

## 先理解数据对象

- **counts**：原始计数矩阵，通常行是基因、列是细胞。零表示这次测序未检出，不一定表示细胞真正完全不表达。
- **data**：通常是标准化后的 log 表达。与 counts 的数值尺度不同，不能随意替换。
- **Seurat 对象**：可以理解为分析资料包，包含表达矩阵、逐细胞信息和降维坐标。
- **meta.data / obj[[]]**：逐细胞信息表，通常每行一个细胞，列保存样本、细胞类型等。
- **barcode / cell_id**：细胞的唯一名字。给细胞添加标签时，必须按名字对齐，不能只看行号。
- **RDS**：R 对象存档，使用 `readRDS()` 恢复；CSV/TSV 是表格，适合用 Excel 查看。

## 最常见的 R 语法

```r
x <- 3                         # 把 3 保存到变量 x
genes <- c("STK31", "HLA-A")    # c() 把多个值组成向量
meta$sample                    # 取 meta 表的 sample 列
mat[genes, , drop = FALSE]      # 取指定基因的行，保留全部细胞列和二维结构
keep <- meta$sample == "P1"     # 为每行得到 TRUE 或 FALSE
meta[keep, ]                   # 普通 data.frame：只保留 TRUE 对应的行
mean(c(TRUE, FALSE, TRUE))      # TRUE=1，FALSE=0，因此结果为 2/3
```

`function(...)` 先定义一段可以复用的操作，调用函数名才会执行函数体。`if` 选择条件分支，`for` 重复处理，`lapply` 对列表逐项运行同一函数。`NA` 表示缺失值；`is.na()` 检查缺失，不能把 NA 当成 0。

`data.table` 的语法与普通表格略有不同：`DT[条件, 计算, by=分组]` 表示筛选并分组计算，`:=` 更新列，`.N` 统计行数。`ggplot` 的 `aes()` 指定“哪个列映射到坐标/颜色/大小”，后面的 `+` 是继续添加图层。

## 本项目特别容易混淆的规则

### STK31-high 的定义随分析版本不同

| 位置 | 实际规则 |
|---|---|
| 02/03/06 的上皮分组函数 | 计算上皮内第 75 百分位；阈值大于零时取 `>= cutoff`，阈值为零时取 `>0` |
| 04 的 UMAP 分组 | 沿用已有 cutoff，但代码使用严格的 `> cutoff`；恰好等于阈值的细胞可能与主流程分组不同 |
| 05 的 STK31-positive | 上皮原始 `counts > 0` |
| 07 的敏感性检查 | 对照 counts>0、全上皮 q75、逐样本上皮 q75，并检查细胞集合是否实际相同 |
| 11/18 的 CNV 恶性上皮分组 | 恶性上皮内 `counts > 0` 为 high，`counts == 0` 为 low |

因此 high 不总是“表达最高的 25%”，不同版本的 high 也不一定是同一群细胞。旧版 `Epithelial` 来自细胞类型注释，不能仅凭名称中的 tumor 就认为已经通过 CNV 确认恶性。

### CNV 的 Normal 是当前工作分组

脚本 11 将 `malignant_high_confidence` 和 `malignant_likely` 归恶性；将 `cnv_normal_like`、`uncertain`、`unresolved_low_quality` 和缺少 CNV 结果的上皮归 `Normal epithelial cells`。后者包含不确定细胞，不表示每一个都已证实正常。原始 CNV 状态与标签来源仍在对象中保留。

样本映射是 **tissue1→P1、tissue2→P2、tissue4→P3、tissue5→P4**。

### 先确认统计单位和比较方向

`avg_log2FC > 0` 通常表示差异检验的第一组比第二组表达更高。调整后 P 值处理了多基因检验问题，但不会自动解决同一患者内细胞相关的问题。

细胞级 Wilcoxon 用细胞作比较；患者配对伪批量先将“患者×分组”的原始计数相加，再用设计矩阵考虑患者效应。脚本 11 的配对分析包括 P1/P2/P4，共六个患者×组别汇总列，并不是六个独立患者。

脚本 07 的 bootstrap 是样本内抽细胞，LOSO 是每次去掉一个样本。二者检查的稳定性不同，均不增加患者数。

### 读图时先看数值尺度

- UMAP：颜色可表示细胞类型或表达，坐标不等于组织中的物理位置。
- 检测比例：表达>0 的细胞数除以当前组总细胞数。
- 伪批量 CPM：组内某基因总计数除以该组所有基因总计数，再乘一百万；与逐细胞标准化后取均值不同。
- z-score 热图：颜色是相对该行均值的偏离，不能直接比较不同基因的绝对表达。
- CellChat：发送→接收的计算推断；未检出不等于真实不存在。人工候选表评分与 CellChat 概率也不同。
- GO：说明候选基因与某过程的关联；上下调混合后富集不能直接解释成过程被激活。
- 虚拟敲除：网络扰动幅度与观察性匹配得到的方向是两种证据，不能替代真实敲除实验。

## 阅读与运行建议

每次只读一个脚本的一个小区段：先问“输入是什么”，再看“筛选了哪些细胞”，最后找 `write.csv`、`fwrite`、`saveRDS` 或 `ggsave` 确认输出去向。先理解函数说明，再进入函数体，最后回到主流程看调用顺序。

本项目 R 分析默认在 Linux 服务器 `/home/zhuweiyu/codex-r` 运行，完整分析依赖外部数据和已安装的包。本地 Python 图/PPT 脚本有 Windows 路径和字体配置。PPT 中部分标题、日期、人数和结论是固定文本，更新数据后需要逐项核对。

这次修改仅补充阅读注释，没有重新运行单细胞分析或更新结果图表。
