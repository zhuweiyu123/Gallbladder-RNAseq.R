#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
out <- file.path(project, "06_feasibility", "stage3")
summary_file <- file.path(out, "stage3_summary.md")

mal <- fread(file.path(out, "malignant_target_gene_detection.csv"))
mal_sample <- fread(file.path(out, "malignant_sample_gene_detection.csv"))
stk_threshold <- fread(file.path(out, "stk31_sample_threshold_summary.csv"))
nk <- fread(file.path(out, "nk_target_gene_detection.csv"))
nk_sample <- fread(file.path(out, "nk_sample_gene_detection.csv"))
kir_threshold <- fread(file.path(out, "kir3dl1_sample_threshold_summary.csv"))
nk_subtype <- fread(file.path(out, "nk_subtype_gene_detection.csv"))
exact <- fread(file.path(out, "primary_matched_stk31_hla_nk.csv"))
structural <- fread(file.path(out, "primary_structural_sensitivity_stage3.csv"))
grid <- fread(file.path(out, "primary_threshold_sensitivity_stage3.csv"))
coverage <- fread(file.path(out, "primary_reference_signal_coverage.csv"))
runtime <- fread(file.path(out, "stage3_runtime_memory.csv"))
validation <- fread(file.path(out, "stage3_validation_checks.csv"))
concentration_mal <- fread(file.path(out, "malignant_signal_concentration.csv"))
concentration_nk <- fread(file.path(out, "nk_signal_concentration.csv"))

fmt_int <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x, digits = 2L) sprintf(paste0("%.", digits, "f%%"), 100 * as.numeric(x))
row_for <- function(data, requested_gene) data[data$gene == requested_gene][1]

mal_table <- paste(
  c(
    "| 基因 | RNA检出细胞/总细胞 | 检出率 | raw UMI总数 | 检出细胞raw count中位数 |",
    "|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(mal)), function(i) sprintf(
      "| %s | %s / %s | %s | %s | %.1f |",
      mal$gene[i], fmt_int(mal$RNA_detected_n[i]), fmt_int(mal$total_cell_n[i]),
      fmt_pct(mal$RNA_detected_fraction[i], 2), fmt_int(mal$total_raw_UMI[i]),
      mal$median_raw_count_among_detected[i]
    ), character(1))
  ), collapse = "\n"
)

nk_table <- paste(
  c(
    "| 基因 | RNA检出细胞/总NK | 检出率 | raw UMI总数 |",
    "|---|---:|---:|---:|",
    vapply(seq_len(nrow(nk)), function(i) sprintf(
      "| %s | %s / %s | %s | %s |",
      nk$gene[i], fmt_int(nk$RNA_detected_n[i]), fmt_int(nk$total_NK_n[i]),
      fmt_pct(nk$RNA_detected_fraction[i], 2), fmt_int(nk$total_raw_UMI[i])
    ), character(1))
  ), collapse = "\n"
)

threshold_table <- paste(
  c(
    "| 集合 | ≥1 | ≥3 | ≥5 | ≥10 | ≥20 |",
    "|---|---:|---:|---:|---:|---:|",
    sprintf(
      "| STK31：全部58个恶性样本 | %s | %s | %s | %s | %s |",
      stk_threshold[scope == "all_malignant_samples" & threshold_RNA_detected_cell_n == 1, sample_n],
      stk_threshold[scope == "all_malignant_samples" & threshold_RNA_detected_cell_n == 3, sample_n],
      stk_threshold[scope == "all_malignant_samples" & threshold_RNA_detected_cell_n == 5, sample_n],
      stk_threshold[scope == "all_malignant_samples" & threshold_RNA_detected_cell_n == 10, sample_n],
      stk_threshold[scope == "all_malignant_samples" & threshold_RNA_detected_cell_n == 20, sample_n]
    ),
    sprintf(
      "| STK31：11个有恶性细胞的原发样本 | %s | %s | %s | %s | %s |",
      stk_threshold[scope == "primary_tumor_malignant_samples" & threshold_RNA_detected_cell_n == 1, sample_n],
      stk_threshold[scope == "primary_tumor_malignant_samples" & threshold_RNA_detected_cell_n == 3, sample_n],
      stk_threshold[scope == "primary_tumor_malignant_samples" & threshold_RNA_detected_cell_n == 5, sample_n],
      stk_threshold[scope == "primary_tumor_malignant_samples" & threshold_RNA_detected_cell_n == 10, sample_n],
      stk_threshold[scope == "primary_tumor_malignant_samples" & threshold_RNA_detected_cell_n == 20, sample_n]
    ),
    sprintf(
      "| KIR3DL1：全部133个NK样本 | %s | %s | %s | %s | %s |",
      kir_threshold[scope == "all_NK_samples" & threshold_RNA_detected_cell_n == 1, sample_n],
      kir_threshold[scope == "all_NK_samples" & threshold_RNA_detected_cell_n == 3, sample_n],
      kir_threshold[scope == "all_NK_samples" & threshold_RNA_detected_cell_n == 5, sample_n],
      kir_threshold[scope == "all_NK_samples" & threshold_RNA_detected_cell_n == 10, sample_n],
      kir_threshold[scope == "all_NK_samples" & threshold_RNA_detected_cell_n == 20, sample_n]
    ),
    sprintf(
      "| KIR3DL1：86个有NK的原发样本 | %s | %s | %s | %s | %s |",
      kir_threshold[scope == "primary_tumor_NK_samples" & threshold_RNA_detected_cell_n == 1, sample_n],
      kir_threshold[scope == "primary_tumor_NK_samples" & threshold_RNA_detected_cell_n == 3, sample_n],
      kir_threshold[scope == "primary_tumor_NK_samples" & threshold_RNA_detected_cell_n == 5, sample_n],
      kir_threshold[scope == "primary_tumor_NK_samples" & threshold_RNA_detected_cell_n == 10, sample_n],
      kir_threshold[scope == "primary_tumor_NK_samples" & threshold_RNA_detected_cell_n == 20, sample_n]
    )
  ), collapse = "\n"
)

kir_subtype <- nk_subtype[gene == "KIR3DL1"][order(-RNA_detected_fraction)]
kir_subtype_table <- paste(
  c(
    "| 作者NK亚群 | 细胞数 | KIR3DL1 RNA检出数 | 检出率 | Wilson 95%区间 |",
    "|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(kir_subtype)), function(i) sprintf(
      "| %s | %s | %s | %s | %s–%s |",
      kir_subtype$nk_subtype[i], fmt_int(kir_subtype$subtype_cell_n[i]),
      fmt_int(kir_subtype$RNA_detected_n[i]), fmt_pct(kir_subtype$RNA_detected_fraction[i], 2),
      fmt_pct(kir_subtype$wilson95_lower[i], 2), fmt_pct(kir_subtype$wilson95_upper[i], 2)
    ), character(1))
  ), collapse = "\n"
)

primary_mal_coverage <- mal_sample[standardized_site == "Primary tumor", .(
  primary_sample_n = .N,
  pseudobulk_nonzero_sample_n = sum(pseudobulk_raw_count > 0),
  pseudobulk_nonzero_fraction = mean(pseudobulk_raw_count > 0)
), by = gene]

receptor_primary <- concentration_nk[
  scope == "primary_tumor_NK_samples" & signal_metric == "pseudobulk_raw_count" &
    gene %in% c("KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRC1")
]
receptor_table <- paste(
  c(
    "| NK受体 | 86个原发NK样本中RNA非零 | 覆盖率 | 全部NK细胞检出率 |",
    "|---|---:|---:|---:|",
    vapply(seq_len(nrow(receptor_primary)), function(i) {
      gene <- receptor_primary$gene[i]
      overall <- row_for(nk, gene)
      sprintf("| %s | %s / 86 | %s | %s |", gene, receptor_primary$nonzero_unit_n[i],
              fmt_pct(receptor_primary$nonzero_unit_fraction[i], 2), fmt_pct(overall$RNA_detected_fraction, 2))
    }, character(1))
  ), collapse = "\n"
)

reference <- exact[structurally_evaluable_reference == TRUE]
reference_members <- paste(reference$sample_name, collapse = "、")
reference_signal <- grid[
  malignant_cell_threshold == 100 & NK_cell_threshold == 100 &
    STK31_detected_cell_threshold == 1 & KIR3DL1_detected_cell_threshold == 1
][1]
strict_signal <- grid[
  malignant_cell_threshold == 100 & NK_cell_threshold == 100 &
    STK31_detected_cell_threshold == 5 & KIR3DL1_detected_cell_threshold == 5
][1]
stk_concentration <- concentration_mal[
  gene == "STK31" & scope == "all_malignant_samples" & signal_metric == "pseudobulk_raw_count"
][1]

runtime_table <- paste(
  c(
    "| 步骤 | 用时 | 峰值RSS | 退出码 |",
    "|---|---:|---:|---:|",
    vapply(seq_len(nrow(runtime)), function(i) sprintf(
      "| %s | %s | %.3f GiB | %s |",
      runtime$component[i], runtime$elapsed[i], runtime$max_rss_gib[i], runtime$exit_status[i]
    ), character(1))
  ), collapse = "\n"
)

lines <- c(
  "# 第三阶段：STK31–HLA-I–NK 基因检出率与样本可行性筛查",
  "",
  "## 执行状态",
  "",
  sprintf("第三阶段已完成并停止在可行性筛查范围内。所有 %d 项跨表与结构校验均通过；未进行相关性、p值、差异表达、score、CellChat、富集分析、重聚类或任何机制推断。", sum(validation$passed)),
  "",
  "本报告中的“RNA检出”统一定义为原始 RNA assay counts > 0。scRNA零值可能来自掉零；RNA检出不等同于蛋白表达或表面受体阳性。",
  "",
  "## 1. 安全提取与数据产品",
  "",
  "- 采用 R 生成并核验manifest，再以 C++/zlib 单次流式扫描全细胞 `matrix.mtx.gz`；没有在R中建立111万细胞对象，也未转为dense matrix。",
  "- 输入核验：32,137基因 × 1,117,245细胞，声明和实际读取均为1,457,063,953个非零值；gzip完整读至EOF，索引均合法，无相邻重复坐标。",
  "- 完整恶性上皮raw-count矩阵：32,137 × 96,942，203,691,975个非零值；保存为标准稀疏Matrix Market三件套 `02_processed_data/malignant_epithelial/10X_counts/`，矩阵压缩后约591 MiB。",
  "- 目标小矩阵：预设8基因加补充HLA-E，共9 × 96,942；保存于 `target_counts/`，并另存约631 KiB的 `malignant_target_raw_counts.rds`。",
  "- 恶性细胞metadata共96,942行，barcode与矩阵列顺序完全一致；`patient_id_source=derived_from_sample_prefix`，不是作者官方患者ID。",
  "- HLA-E是为KLRC1/NKG2A可行性补充核查，不改变预设8基因主面板，也未构建任何综合score。",
  "",
  "## 2. 恶性上皮目标基因检出",
  "",
  mal_table,
  "",
  sprintf("STK31在96,942个恶性细胞中检出2,885个（2.98%%），总raw UMI为3,293；检出细胞中位数为1。这说明它在单细胞层面稀疏，不适合据此定义STK31-high/low，但仍可检查样本级raw pseudobulk。58个恶性样本中31个pseudobulk非零（53.45%%），11个原发恶性样本中7个非零（63.64%%）。全部恶性样本的STK31 raw counts有%s由最高一个样本贡献、%s由最高三个样本贡献，存在明显样本集中，但未达到预设top-1 >50%%“高度集中”警戒。", fmt_pct(stk_concentration$top1_share, 2), fmt_pct(stk_concentration$top3_share, 2)),
  "",
  "HLA-A/B/C和B2M的细胞级检出率分别为76.64%、73.96%、75.07%和93.18%，明显高于STK31；NLRC5、TAP1、TAP2较稀疏，但在11个原发恶性样本中的pseudobulk非零覆盖分别为9/11、10/11、9/11。因此，经典HLA-I表达/抗原呈递基因适合进入样本级pseudobulk描述与后续建模，不能据此证明特定等位基因配体关系。",
  "",
  "## 3. STK31与KIR3DL1样本覆盖",
  "",
  threshold_table,
  "",
  "这些阈值结果只描述可测信号覆盖，不作为正式分析的纳入/排除规则，避免按暴露或受体信号预筛样本。",
  "",
  "## 4. NK目标基因检出",
  "",
  nk_table,
  "",
  "KIR3DL1在28,266个NK中检出493个（1.74%），总raw UMI 556，中位raw count为1，属于稀疏受体。其pseudobulk在133个NK样本中80个非零（60.15%），在86个原发NK样本中53个非零（61.63%）。它可以保留为候选辅助轴和样本级稀疏计数变量，但不适合作为单细胞主分组或唯一机制支柱。KIR高度同源且多态，3′ scRNA比对与掉零可能低估KIR3DL1；不能将零值解释为NK不表达。",
  "",
  receptor_table,
  "",
  "KLRC1最稳定；KIR2DL3也比KIR3DL1广泛；KIR2DL1与KIR3DL1样本覆盖接近但细胞级略高。广义HLA-I抗原呈递与NK抑制性受体共存的测量可行性优于单独HLA-B–KIR3DL1轴。",
  "",
  "### KIR3DL1按作者NK亚群的观察到的RNA检出比例",
  "",
  kir_subtype_table,
  "",
  "这里只做比例排序和Wilson区间，没有富集检验，也未改变作者亚群注释。最高观察比例为NK_C10_HLA-DRA（58/691，8.39%），其次为NK_C0_FCGR3A（229/7,055，3.25%）和NK_C5_ISG15（55/1,760，3.12%）；小亚群的不确定性以区间展示。",
  "",
  "## 5. 原发肿瘤精确配对与阈值敏感性",
  "",
  "- sample_info中有86个scRNA原发肿瘤样本；恶性端和NK端只通过完全相同的sample_name做inner match，得到11个精确同一样本配对、11个派生患者。没有将同患者不同样本拼成配对。",
  sprintf("- 预先锁定的结构参考条件为：同一原发样本、恶性细胞≥100、NK细胞≥100、两端raw library>0。仅%d个独立样本/患者满足：%s。", nrow(reference), reference_members),
  sprintf("- 在这4个结构可评估样本中，STK31 pseudobulk为4/4非零，KIR3DL1为3/4非零；若另外要求STK31检出细胞≥1且KIR3DL1检出NK≥1，只剩%d个：%s。", reference_signal$eligible_sample_n, reference_signal$sample_names),
  sprintf("- 若在同一结构条件下要求STK31检出≥5且KIR3DL1检出≥5，只剩%d个：%s。", strict_signal$eligible_sample_n, ifelse(nchar(strict_signal$sample_names), strict_signal$sample_names, "无")),
  "- 仅按细胞数的16组合中，可评估独立患者为3–7个；完整4×4×4×4=256个组合均已保留（包括0结果），其signal-supported患者数为0–4。",
  "- 由于参考结构集合只有4个患者，患者级STK31–HLA-I–NK关联不适合在此数据中作为稳健/确证分析；至多作为非常小样本的描述性或pilot结果。不能通过挑选宽松或有利阈值来扩大结论。",
  "",
  "## 6. 八个必须回答的问题",
  "",
  "1. **STK31是否足够？** 单细胞层面不够稳定（2.98%），不支持STK31-high/low；样本级raw pseudobulk可测（全部恶性样本31/58非零，原发恶性样本7/11非零），可作为探索性计数变量，但样本集中且原发精确配对太少。",
  "2. **HLA-B及广义HLA-I是否足够？** 足够用于表达层面的样本级pseudobulk。HLA-A/B/C、B2M细胞级检出高，原发恶性样本均11/11非零；NLRC5/TAP1/TAP2也有较高样本覆盖。",
  "3. **KIR3DL1是否足够？** 单细胞层面稀疏（493/28,266，1.74%）；样本级有一定覆盖（80/133、原发53/86），但精确配对参考集合仅3/4非零。",
  "4. **KIR3DL1定位？** 只能作为候选辅助轴/敏感性指标，不宜作为单细胞重点轴或唯一机制轴。",
  "5. **其他抑制受体是否更稳定？** KLRC1明确更稳定；KIR2DL3也更广；KIR2DL1仅略高于KIR3DL1。KIR2DL2因feature不存在，已完全排除，未填0。",
  "6. **原发肿瘤最终有多少独立样本同时具备结构与信号？** 结构参考条件下4个；若再描述性要求STK31≥1且KIR3DL1≥1，则3个；若两者均≥5，则1个。后两类不是正式纳入规则。",
  "7. **KIR3DL1不足时，广义轴是否可行？** 测量层面可行：KLRC1在85/86原发NK样本非零，KIR2DL3为72/86，且补充HLA-E在原发恶性样本中11/11非零。但这只能支持广义HLA-I抗原呈递与NK抑制受体共存的探索，不能证明特异受体–配体机制。",
  "8. **是否值得进入第四阶段？** **有条件值得。** 适合把完整恶性稀疏矩阵用于STK31/HLA-I的样本级pseudobulk和描述性分析，把NK抑制受体作为辅助证据；不适合宣称公共数据足以验证完整STK31→HLA-B→KIR3DL1患者级机制。第四阶段若进行，应预先降级为探索性/pilot，并把更大或独立队列、HLA分型和蛋白/功能实验作为关键验证。",
  "",
  "## 7. 科学解释边界",
  "",
  "- 总HLA-B RNA不能替代HLA-Bw4分型；KIR3DL1的直接配体特异性无法由当前数据确认。部分HLA-A也可能携带Bw4，但同样需要分型。",
  "- 总HLA-C RNA不能区分HLA-C1/C2，因而不能证明KIR2DL3/KIR2DL1配体关系。",
  "- KLRC1/NKG2A的相关配体是HLA-E；本阶段仅补充核查HLA-E feature与表达。KLRD1/CD94单独不等于抑制性信号。",
  "- RNA检出不等于蛋白、受体表面量或功能抑制；没有计算相关性、score或机制效应。",
  "- patient_id来自样本名前缀推导，仍不是作者明确提供的官方患者ID；在正式发表前需由原作者补充材料或数据说明进一步确认。",
  "",
  "## 8. 性能、资源与验证",
  "",
  runtime_table,
  "",
  "服务器当前约503 GiB内存、约533 GiB可用磁盘。Stage 3最高峰值RSS约1.915 GiB（读取NK.RDS）；14.57亿非零值流式提取仅约0.152 GiB。完整恶性稀疏矩阵约591 MiB，因此当前服务器资源充足，但仍不建议把完整全细胞矩阵读入R或转成dense。",
  "",
  sprintf("验证结果：%d/%d项通过，0项失败；同时对两个输出matrix.mtx.gz执行了`gzip -t`。NK第一次运行因脚本中标量/向量`fifelse`长度错误在输出前停止，修复后重跑成功；失败日志和time文件保留，最终结果只来自成功运行。", sum(validation$passed), nrow(validation)),
  "",
  "## 9. 主要输出",
  "",
  "- 规定的11个核心输出均已生成，包括完整256行阈值网格和本报告。",
  sprintf("- 另生成结构阈值、信号集中度、全86个原发样本区室可用性、HLA-E补充核查、运行资源和%d项验证明细。", nrow(validation)),
  "- 所有R/bash/C++脚本保存在`04_scripts/`；预先锁定的决策规则保存在`00_docs/stage3_decision_rules.md`。",
  "",
  "## 10. 尚待核实",
  "",
  "- **待核实：** 样本名前缀是否完全等同于作者定义的患者ID；本阶段仅使用明确标记的派生ID。",
  "- **待核实：** 各患者HLA-Bw4、HLA-C1/C2及其他HLA等位基因背景；当前文件未提供可直接使用的HLA分型。",
  "- **待核实：** KIR3DL1/KIR2DL1/KIR2DL3的蛋白表面表达、克隆性和功能抑制状态；scRNA counts不能替代。",
  "- **待核实：** 独立队列能否提供足够患者级样本量验证STK31–HLA-I–NK关联。",
  "",
  "---",
  "",
  "第三阶段到此停止；未自动进入第四阶段。"
)

writeLines(lines, summary_file, useBytes = TRUE)
cat(paste0("WROTE ", summary_file, "\n"))
