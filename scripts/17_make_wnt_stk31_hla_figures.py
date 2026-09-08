#!/usr/bin/env python3
# ========================================================================
# 【中文阅读指南】用已有 CSV 制作 WNT 与 STK31/HLA 专题图
# 输入：脚本 16 的 WNT 汇总、脚本 07 的样本候选轴表、脚本 11 的恶性上皮高低差异表。
# 流程：pandas 读表和筛选 → Pillow 在画布上绘制 → 保存 PNG/PDF 及对应绘图数据。
# 输出：results/stk31_wnt_hla_focused_figures；字体使用 Windows 的 Arial 路径。
# 图 2 的样本面板来自旧 refined 上皮定义，基因面板来自 CNV 恶性上皮定义，两者数据范围不同。
# 这份脚本负责展示既有统计结果；标题和样本展示顺序是固定文本，换数据时需重新核对。
# Python 入门：def 定义函数；缩进表示代码归属；字典保存键值对；Path 管理文件路径。
# 先读顶部输入路径与函数说明，再看文件末尾入口；不要把图中文字当作自动生成的统计结论。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
"""Create two focused, publication-ready figures from existing analysis outputs.

Figure 1: WNT CellChat source-target probabilities in datasets with detected WNT.
Figure 2: Association between STK31-high status and HLA-I expression.

The script does not rerun CellChat or differential expression. It only visualizes
previously exported tables and writes the exact plotting data beside the figures.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
# 【图像工具包】Pillow 提供画布、文字和线条操作；这里按像素排版，最终保存成图片。
from PIL import Image, ImageDraw, ImageFont


# 【项目根目录】后面的输入输出路径以此为起点；相对脚本定位与写死本机路径的可移植性不同。
PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "stk31_wnt_hla_focused_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WNT_INPUT = (
    PROJECT_DIR
    / "results"
    / "qc4_local_sample_cellchat_pathways"
    / "00_combined"
    / "WNT_source_target_probabilities.csv"
)
HLA_SAMPLE_INPUT = (
    PROJECT_DIR
    / "results"
    / "refined_sample_level_validation_v2"
    / "sample_level_candidate_axis_summary.csv"
)
HLA_GENE_INPUT = (
    PROJECT_DIR
    / "results"
    / "merged_cnv_binary_annotation"
    / "differential_expression"
    / "markers_stk31_high_vs_low_malignant_epithelial.csv"
)

FONT_REGULAR = Path(r"C:\Windows\Fonts\arial.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")


# 【函数：font】按 size 和 bold 读取 Arial 常规/粗体字体；返回 Pillow 绘字需要的字体对象。
# 字体路径在文件顶部配置，文件不存在会在加载时报错。
def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REGULAR), size=size)


# 【函数：save_png_pdf】将 Pillow 画布保存为 PNG，再转 RGB 保存 PDF；输出分辨率设为 320 dpi。
# 这个 PDF 包含栅格图片，放大时不会变成可无限缩放的矢量图。
def save_png_pdf(image: Image.Image, stem: Path) -> None:
    image.save(stem.with_suffix(".png"), dpi=(320, 320), optimize=True)
    image.convert("RGB").save(stem.with_suffix(".pdf"), "PDF", resolution=320)


# 【函数：draw_centered】先测量文字边界，再偏移起点，让文字在给定 xy 附近居中。
def draw_centered(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str,
                  text_font: ImageFont.FreeTypeFont, fill: str = "#111827") -> None:
    box = draw.textbbox((0, 0), text, font=text_font)
    draw.text((xy[0] - (box[2] - box[0]) / 2, xy[1] - (box[3] - box[1]) / 2),
              text, font=text_font, fill=fill)


# 【函数：draw_rotated_text】先在透明小画布上写文字，再旋转并贴到主图，用于竖排轴标签。
def draw_rotated_text(image: Image.Image, center: tuple[int, int], text: str,
                      text_font: ImageFont.FreeTypeFont, angle: float,
                      fill: str = "#374151") -> None:
    box = text_font.getbbox(text)
    width = box[2] - box[0] + 30
    height = box[3] - box[1] + 30
    # 【像素画布】width/height 指像素；RGBA 有透明通道，draw 负责在画布上添加元素。
    layer = Image.new("RGBA", (width, height), (255, 255, 255, 0))
    layer_draw = ImageDraw.Draw(layer)
    layer_draw.text((15 - box[0], 15 - box[1]), text, font=text_font, fill=fill)
    rotated = layer.rotate(angle, expand=True, resample=Image.Resampling.BICUBIC)
    image.alpha_composite(rotated, (int(center[0] - rotated.width / 2),
                                    int(center[1] - rotated.height / 2)))


# 【函数：interpolate_color】按 value 在相邻颜色节点之间插值，生成连续颜色梯度。
# 颜色只是数值编码，具体编码哪个指标以调用处为准。
def interpolate_color(stops: list[tuple[float, str]], value: float) -> str:
    value = max(0.0, min(1.0, value))
    for idx in range(len(stops) - 1):
        left, right = stops[idx], stops[idx + 1]
        if left[0] <= value <= right[0]:
            frac = (value - left[0]) / (right[0] - left[0])
            lrgb = tuple(int(left[1][i:i + 2], 16) for i in (1, 3, 5))
            rrgb = tuple(int(right[1][i:i + 2], 16) for i in (1, 3, 5))
            rgb = tuple(round(a + frac * (b - a)) for a, b in zip(lrgb, rrgb))
            return "#" + "".join(f"{channel:02X}" for channel in rgb)
    return stops[-1][1]


# 【函数：make_wnt_figure】读取 WNT 汇总表，按预设样本和细胞类型筛选发送/接收组合，再绘制网络图。
# 先保存实际绘图 CSV，便于从图追溯数据；样本候选顺序在函数中明确列出。
def make_wnt_figure() -> None:
    # 【读表】pandas 将 CSV 转成 DataFrame；后续用列名选择指标，用布尔条件筛选需要的行。
    data = pd.read_csv(WNT_INPUT)
    # 【数值转换】无法解析的文本在 errors="coerce" 下变成 NaN；后续 fillna(0) 会把缺失按零处理。
    data["probability"] = pd.to_numeric(data["probability"], errors="coerce").fillna(0.0)
    detected = (
        data.groupby("dataset", as_index=False)["probability"]
        .sum()
        .query("probability > 0")["dataset"]
        .tolist()
    )
    # 【固定展示范围】仅从这里列出的 tissue1/2/5 中选已检出 WNT 的数据集，其他样本不会自动加入。
    dataset_order = [x for x in ["tissue1", "tissue2", "tissue5"] if x in detected]
    group_order = [
        "Malignant epithelial cells",
        "Normal epithelial cells",
        "Fibroblast/stromal cells",
        "Endothelial cells",
    ]
    plot_data = data[
        data["dataset"].isin(dataset_order)
        & data["source"].isin(group_order)
        & data["target"].isin(group_order)
    ].copy()
    # 【绘图数据留档】导出实际用于图中的子表；index=False 避免额外写入 DataFrame 行索引。
    plot_data.to_csv(OUT_DIR / "Figure1_WNT_plot_data.csv", index=False)

    width, height = 3600, 1850
    image = Image.new("RGBA", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((150, 65), "WNT signaling source-target networks",
              font=font(82, bold=True), fill="#111827")
    draw.text(
        (150, 165),
        "CellChat probabilities in samples with detected WNT signaling",
        font=font(42),
        fill="#4B5563",
    )
    draw.text(
        (150, 225),
        "WNT was not detected in QC4_public or tissue4 under the current inference settings.",
        font=font(34),
        fill="#6B7280",
    )

    abbreviated = {
        "Malignant epithelial cells": "Malignant epi.",
        "Normal epithelial cells": "Normal epi.",
        "Fibroblast/stromal cells": "Fibroblast/stromal",
        "Endothelial cells": "Endothelial",
    }
    heat_top = 390
    cell = 185
    panel_x = [470, 1360, 2250]
    max_probability = float(plot_data["probability"].max())
    color_stops = [
        (0.00, "#F7FBFF"),
        (0.30, "#C6DBEF"),
        (0.60, "#6BAED6"),
        (1.00, "#08519C"),
    ]

    for panel_index, dataset in enumerate(dataset_order):
        x0 = panel_x[panel_index]
        draw_centered(draw, (x0 + 2 * cell, 335), dataset,
                      font(48, bold=True), fill="#1F2937")
        subset = plot_data[plot_data["dataset"] == dataset]
        lookup = {(row.source, row.target): float(row.probability)
                  for row in subset.itertuples(index=False)}
        for row_idx, target in enumerate(group_order):
            for col_idx, source in enumerate(group_order):
                value = lookup.get((source, target), 0.0)
                scaled = math.sqrt(value / max_probability) if max_probability > 0 else 0.0
                fill = interpolate_color(color_stops, scaled)
                left = x0 + col_idx * cell
                top = heat_top + row_idx * cell
                draw.rectangle((left, top, left + cell, top + cell), fill=fill,
                               outline="#D1D5DB", width=3)
                if value > 0:
                    label = f"{value * 1000:.2f}"
                    text_fill = "white" if scaled > 0.62 else "#111827"
                    draw_centered(draw, (left + cell / 2, top + cell / 2), label,
                                  font(28, bold=scaled > 0.62), fill=text_fill)
        if panel_index == 0:
            for row_idx, target in enumerate(group_order):
                y = heat_top + (row_idx + 0.5) * cell
                label = abbreviated[target]
                box = draw.textbbox((0, 0), label, font=font(31))
                draw.text((x0 - 35 - (box[2] - box[0]), y - 18), label,
                          font=font(31), fill="#374151")
        for col_idx, source in enumerate(group_order):
            x = x0 + (col_idx + 0.5) * cell
            draw_rotated_text(image, (int(x), heat_top + 4 * cell + 135),
                              abbreviated[source], font(30), 42)

    draw_rotated_text(image, (105, heat_top + 2 * cell), "WNT receiver",
                      font(40, bold=True), 90, fill="#1F2937")
    draw_centered(draw, (1790, 1485), "WNT sender", font(42, bold=True), fill="#1F2937")
    draw_centered(draw, (1790, 1560), "Cell values: probability x 10^3",
                  font(31), fill="#6B7280")

    legend_x, legend_y, legend_w, legend_h = 3210, 480, 75, 650
    for step in range(legend_h):
        frac = 1 - step / (legend_h - 1)
        fill = interpolate_color(color_stops, frac)
        draw.line((legend_x, legend_y + step, legend_x + legend_w, legend_y + step),
                  fill=fill, width=2)
    draw.rectangle((legend_x, legend_y, legend_x + legend_w, legend_y + legend_h),
                   outline="#9CA3AF", width=2)
    draw.text((3120, 380), "CellChat probability", font=font(33, bold=True), fill="#1F2937")
    draw.text((3150, 425), "sqrt color scale", font=font(27), fill="#6B7280")
    for tick in [0.0, 0.005, 0.010, 0.015, max_probability]:
        frac = math.sqrt(tick / max_probability) if max_probability > 0 else 0.0
        y = legend_y + legend_h * (1 - frac)
        draw.line((legend_x + legend_w, y, legend_x + legend_w + 14, y),
                  fill="#374151", width=2)
        draw.text((legend_x + legend_w + 24, y - 17), f"{tick:.3f}",
                  font=font(27), fill="#374151")

    draw.text(
        (150, 1750),
        "Source: existing CellChat output; identical inference settings across datasets.",
        font=font(29),
        fill="#6B7280",
    )
    save_png_pdf(image, OUT_DIR / "Figure1_WNT_source_target_network")


# 【函数：make_hla_figure】分别读取样本候选轴表和恶性上皮差异表，绘制 HLA 的样本与基因面板。
# 两种输入来自不同分析定义；先保存各自绘图数据再排版，避免混淆统计单位。
def make_hla_figure() -> None:
    # 【样本面板来源】这里读取脚本 07 的 refined 上皮候选轴摘要；与后面的 CNV 恶性上皮差异表不是同一细胞集合。
    sample_data = pd.read_csv(HLA_SAMPLE_INPUT)
    sample_data = sample_data[sample_data["axis"] == "MHC-I/HLA_epi_to_NK"].copy()
    sample_columns = [
        "sample", "n_high", "n_low",
        "sender_ligand_mean_for_high_context",
        "sender_ligand_mean_for_low_context",
        "delta_interaction_high_minus_low",
    ]
    sample_data = sample_data[sample_columns]
    sample_data.to_csv(OUT_DIR / "Figure2_STK31_HLA_sample_plot_data.csv", index=False)

    # 【基因面板来源】这里读取脚本 11 的 CNV 恶性上皮 high/low 差异结果，保留原有 P 值和 fold change。
    gene_data = pd.read_csv(HLA_GENE_INPUT)
    hla_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F"]
    gene_data = gene_data[gene_data["gene"].isin(hla_genes)].copy()
    gene_data["gene"] = pd.Categorical(gene_data["gene"], categories=hla_genes, ordered=True)
    gene_data = gene_data.sort_values("gene")
    gene_data.to_csv(OUT_DIR / "Figure2_STK31_HLA_gene_plot_data.csv", index=False)

    width, height = 3600, 1900
    image = Image.new("RGBA", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((150, 60), "STK31-high epithelial cells show lower HLA-I expression trends",
              font=font(76, bold=True), fill="#111827")
    draw.text(
        (150, 160),
        "Within-sample composite expression and pooled gene-level effects; association, not causal proof",
        font=font(40),
        fill="#4B5563",
    )

    # Panel A: paired within-sample HLA-I composite expression.
    left_x0, left_x1 = 520, 1690
    plot_top, plot_bottom = 550, 1370
    draw.text((150, 290), "A", font=font(62, bold=True), fill="#111827")
    draw.text((235, 295), "HLA-I/B2M expression by STK31 group",
              font=font(46, bold=True), fill="#111827")
    draw.text((235, 355), "Paired within each tissue sample", font=font(31), fill="#6B7280")
    x_min, x_max = 1.25, 2.05
    ticks = [1.3, 1.5, 1.7, 1.9, 2.0]
    for tick in ticks:
        x = left_x0 + (tick - x_min) / (x_max - x_min) * (left_x1 - left_x0)
        draw.line((x, plot_top, x, plot_bottom), fill="#E5E7EB", width=3)
        draw_centered(draw, (x, plot_bottom + 55), f"{tick:.1f}", font(30), fill="#4B5563")
    sample_order = ["tissue1", "tissue2", "tissue4", "tissue5"]
    sample_lookup = {row.sample: row for row in sample_data.itertuples(index=False)}
    y_positions = [650, 835, 1020, 1205]
    high_color, low_color = "#D55E00", "#0072B2"
    lower_count = 0
    for sample, y in zip(sample_order, y_positions):
        row = sample_lookup[sample]
        high = float(row.sender_ligand_mean_for_high_context)
        low = float(row.sender_ligand_mean_for_low_context)
        lower_count += int(high < low)
        high_x = left_x0 + (high - x_min) / (x_max - x_min) * (left_x1 - left_x0)
        low_x = left_x0 + (low - x_min) / (x_max - x_min) * (left_x1 - left_x0)
        draw.text((235, y - 23), sample, font=font(34, bold=True), fill="#374151")
        draw.line((high_x, y, low_x, y), fill="#9CA3AF", width=12)
        draw.ellipse((low_x - 22, y - 22, low_x + 22, y + 22), fill=low_color,
                     outline="white", width=4)
        draw.ellipse((high_x - 22, y - 22, high_x + 22, y + 22), fill=high_color,
                     outline="white", width=4)
        delta = high - low
        draw.text((1720, y - 23), f"Delta={delta:+.3f}", font=font(29), fill="#4B5563")
    draw.line((left_x0, plot_bottom, left_x1, plot_bottom), fill="#4B5563", width=3)
    draw_centered(draw, ((left_x0 + left_x1) / 2, plot_bottom + 120),
                  "Mean normalized HLA-A/B/C/E + B2M expression",
                  font(32, bold=True), fill="#374151")
    draw.ellipse((490, 440, 525, 475), fill=high_color)
    draw.text((540, 437), "STK31-high", font=font(30), fill="#374151")
    draw.ellipse((820, 440, 855, 475), fill=low_color)
    draw.text((870, 437), "STK31-low", font=font(30), fill="#374151")
    draw.text((235, 1535), f"{lower_count}/4 tissues: lower HLA-I/B2M expression in STK31-high cells",
              font=font(34, bold=True), fill="#1F2937")

    # Divider.
    draw.line((1900, 300, 1900, 1640), fill="#D1D5DB", width=4)

    # Panel B: pooled gene-level effects in latest CNV-defined malignant cells.
    right_x0, right_x1 = 2300, 3400
    draw.text((1980, 290), "B", font=font(62, bold=True), fill="#111827")
    draw.text((2065, 295), "HLA gene effects in CNV-defined malignant cells",
              font=font(42, bold=True), fill="#111827")
    draw.text((2065, 355), "STK31-high versus STK31-low", font=font(31), fill="#6B7280")
    bx_min, bx_max = -0.85, 0.10
    bticks = [-0.8, -0.6, -0.4, -0.2, 0.0]
    btop, bbottom = 550, 1370
    for tick in bticks:
        x = right_x0 + (tick - bx_min) / (bx_max - bx_min) * (right_x1 - right_x0)
        draw.line((x, btop, x, bbottom), fill="#E5E7EB" if tick != 0 else "#6B7280",
                  width=3 if tick != 0 else 5)
        draw_centered(draw, (x, bbottom + 55), f"{tick:.1f}", font(30), fill="#4B5563")
    gene_y = [630, 790, 950, 1110, 1270]
    zero_x = right_x0 + (0 - bx_min) / (bx_max - bx_min) * (right_x1 - right_x0)
    all_negative = True
    for row, y in zip(gene_data.itertuples(index=False), gene_y):
        effect = float(row.avg_log2FC)
        all_negative = all_negative and effect < 0
        effect_x = right_x0 + (effect - bx_min) / (bx_max - bx_min) * (right_x1 - right_x0)
        draw.text((2065, y - 25), str(row.gene), font=font(36, bold=True), fill="#374151")
        draw.line((effect_x, y, zero_x, y), fill="#3B82F6", width=26)
        draw.ellipse((effect_x - 24, y - 24, effect_x + 24, y + 24), fill="#1D4ED8",
                     outline="white", width=4)
        draw.text((effect_x - 125, y - 22), f"{effect:+.2f}", font=font(29), fill="#1E3A8A")
    draw.line((right_x0, bbottom, right_x1, bbottom), fill="#4B5563", width=3)
    draw_centered(draw, ((right_x0 + right_x1) / 2, bbottom + 120),
                  "Average log2 fold change", font(32, bold=True), fill="#374151")
    draw_centered(draw, ((right_x0 + right_x1) / 2, bbottom + 165),
                  "negative = lower in STK31-high", font(29), fill="#6B7280")
    if all_negative:
        draw.text((2065, 1535), "5/5 HLA genes: negative pooled effect",
                  font=font(34, bold=True), fill="#1F2937")
    draw.text((2065, 1590), "All FDR-adjusted P values = 1.0 (trend only)",
              font=font(30), fill="#B45309")

    draw.text(
        (150, 1770),
        "Panel A: refined sample-level analysis. Panel B: latest CNV-binary malignant epithelial differential expression.",
        font=font(28),
        fill="#6B7280",
    )
    save_png_pdf(image, OUT_DIR / "Figure2_STK31_HLA_relationship")


# 【函数：write_readme】将专题图的数据来源与解读说明写到结果目录，和图片一起保存。
def write_readme() -> None:
    text = """Focused WNT and STK31-HLA figures
==================================

Figure1_WNT_source_target_network
- CellChat WNT source-target probabilities for tissue1, tissue2 and tissue5.
- QC4_public and tissue4 had no WNT detected under the current settings.

Figure2_STK31_HLA_relationship
- Panel A: within-sample HLA-I/B2M composite expression in STK31-high versus
  STK31-low epithelial cells from the refined sample-level validation.
- Panel B: pooled HLA gene log2 fold changes in the latest CNV-defined malignant
  epithelial comparison.
- This figure supports an association/trend. It does not establish causality,
  and the pooled HLA gene effects are not FDR-significant.

Each figure is supplied as PNG and PDF. Exact plot data are supplied as CSV.
"""
    (OUT_DIR / "README.txt").write_text(text, encoding="utf-8")


# 【脚本入口】直接运行本文件时执行下面的主函数；作为模块导入时不自动执行这段入口。
if __name__ == "__main__":
    make_wnt_figure()
    make_hla_figure()
    write_readme()
    print(f"Wrote focused figures to: {OUT_DIR}")
