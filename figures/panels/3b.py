import matplotlib
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('3B')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '3B'

d = pd.read_csv(HERE / 'Figure_3B_transcript_data.csv')
g = pd.read_csv(HERE / 'Figure_3B_parent_gene_data.csv').iloc[0]
ink = '#202938'
gray = '#63707D'
orange = '#D55E00'
blue = '#0072B2'
utr = 'white'
cds = '#52606D'
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 9,
        'text.color': ink,
        'axes.labelcolor': ink,
        'axes.edgecolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.75,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
    }
)
fig = plt.figure(figsize = (5, 6), dpi = 150)
fig.text(
    0.035, 0.963, 'B', ha = 'left', va = 'center', fontsize = 14, fontweight = 'bold'
)
title = fig.text(
    0.14,
    0.963,
    r'Transcript analysis identifies higher P/M for $\it{DLK1}$-201 in ASD',
    ha = 'left',
    va = 'center',
    fontsize = 11,
)
left = fig.add_axes([0.18, 0.235, 0.275, 0.60])
right = fig.add_axes([0.655, 0.235, 0.205, 0.60])
share = fig.add_axes([0.48, 0.235, 0.135, 0.60])
share.set_xlim(0, 80)
share.set_ylim(4.60, -0.60)
share.axis('off')
share.set_title('Input share\n(%)', fontsize = 8, pad = 4)
left.set_title('Transcript structure', fontsize = 9, pad = 13)
right.set_title('P/M change', fontsize = 9, pad = 13)
fig.text(
    0.975, 0.865, r'$P_{\mathrm{BH}}$', ha = 'right', va = 'center', fontsize = 9
)
for i, r in d.iterrows():
    share.barh(
        i, r.inputSharePercent, height = 0.12, color = '#9AA6B0', edgecolor = 'none'
    )
    share.text(
        43,
        i,
        f'{r.inputSharePercent:.1f}',
        ha = 'left',
        va = 'center',
        fontsize = 7.5,
    )
    start = 0
    for col, color in [
        ('fivePrimeUtrLength', utr),
        ('codingSequenceLength', cds),
        ('threePrimeUtrLength', utr),
    ]:
        length = float(r[col])
        left.barh(
            i,
            length / 1000,
            left = start / 1000,
            height = 0.19 if col == 'codingSequenceLength' else 0.07,
            color = color,
            edgecolor = 'black',
            linewidth = 0.75,
        )
        start += length
    color = orange if r.log2FoldChange > 0 else blue
    right.errorbar(
        r.log2FoldChange,
        i,
        xerr = 1.95996398454 * r.log2FoldChangeSE,
        fmt = 'o',
        color = color,
        ecolor = gray,
        markerfacecolor = color if r.adjustedPValue < 0.05 else 'white',
        markeredgecolor = 'black' if r.adjustedPValue < 0.05 else color,
        markeredgewidth = 0.8,
        markersize = 5.2,
        capsize = 2,
        lw = 0.8,
    )
    fig.text(
        0.975,
        0.235 + 0.60 * (4.60 - i) / 5.20,
        f'{r.adjustedPValue:.3g}',
        ha = 'right',
        va = 'center',
        fontsize = 8,
    )
left.set_yticks(
    range(5), [r.transcriptName for r in d.itertuples()], fontsize = 8.5
)
left.tick_params(axis = 'y', length = 0, pad = 5)
left.spines['left'].set_visible(False)
left.set_xlim(0, 4.9)
left.set_xticks([0, 2, 4])
left.set_xlabel('Spliced length (kb)', labelpad = 7)
for ax in [left, right]:
    ax.set_ylim(4.60, -0.60)
    ax.tick_params(axis = 'x', labelsize = 8)
    ax.set_axisbelow(True)
    for yi in range(5):
        ax.axhline(yi, color = '#E7EBEF', lw = 0.55, zorder = 0)
right.set_xlim(-5.5, 5)
right.set_xticks([-4, 0, 4])
right.set_yticks([])
right.spines['left'].set_visible(False)
right.axvline(0, color = gray, lw = 0.7, zorder = 0)
right.set_xlabel(r'$\log_2$(ASD/Unaffected)', labelpad = 7)
legend1 = fig.legend(
    handles = [
        Patch(facecolor = utr, edgecolor = 'black', linewidth = 0.4, label = 'UTR'),
        Patch(facecolor = cds, edgecolor = 'black', linewidth = 0.4, label = 'CDS'),
    ],
    loc = 'lower center',
    bbox_to_anchor = (0.3175, 0.108),
    ncol = 2,
    frameon = False,
    fontsize = 8,
    handlelength = 1.2,
    columnspacing = 0.8,
    borderpad = 0.2,
)
legend2 = fig.legend(
    handles = [
        Line2D([], [], color = orange, lw = 1.4, label = 'Higher P/M in ASD'),
        Line2D([], [], color = blue, lw = 1.4, label = 'Lower P/M in ASD'),
        Line2D(
            [],
            [],
            marker = 'o',
            color = ink,
            linestyle = 'none',
            markeredgecolor = 'black',
            markersize = 5,
            label = r'$P_{\mathrm{BH}}<0.05$',
        ),
        Line2D(
            [],
            [],
            marker = 'o',
            color = ink,
            markerfacecolor = 'white',
            linestyle = 'none',
            markeredgecolor = 'black',
            markersize = 5,
            label = r'$P_{\mathrm{BH}}\geq0.05$',
        ),
    ],
    loc = 'lower center',
    bbox_to_anchor = (0.5, 0.012),
    ncol = 2,
    frameon = False,
    fontsize = 8.3,
    handlelength = 1.3,
    columnspacing = 1.8,
    handletextpad = 0.5,
    labelspacing = 0.6,
    borderpad = 0.2,
)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
while title.get_window_extent(renderer).x1 > fig.bbox.width * 0.975:
    title.set_fontsize(title.get_fontsize() - 0.2)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
outside = []
for t in fig.findobj(Text):
    if t.get_visible() and t.get_text().strip():
        b = t.get_window_extent(renderer)
        if (
            b.x0 < 0
            or b.y0 < 0
            or b.x1 > fig.bbox.width
            or b.y1 > fig.bbox.height
        ):
            outside.append(t.get_text())
assert not outside, outside
assert not legend1.get_window_extent(renderer).overlaps(
    legend2.get_window_extent(renderer)
)
assert (
    abs(
        (
            legend1.get_window_extent(renderer).x0
            + legend1.get_window_extent(renderer).x1
        )
        / 2
        - (left.bbox.x0 + left.bbox.x1) / 2
    )
    < 1
)
stem = 'Figure_3B_DLK1_transcript_effects_input_share_5x6'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 3600)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((850, 1020))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 432))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 432)
c.showPage()
c.save()
