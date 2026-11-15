import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2A')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2A'

ROOT = HERE.parent
source = PROJECT_ROOT / 'results/inputDifferentialExpression.tsv'
before = hashlib.sha256(source.read_bytes()).hexdigest()
df = pd.read_csv(source, sep = '\t')
sig = df.adjustedPValue.lt(0.05)
assert len(df) == 21556 and sig.sum() == 305
assert (sig & df.log2FoldChange.gt(0)).sum() == 135
assert (sig & df.log2FoldChange.lt(0)).sum() == 170
ink = '#202938'
muted = '#63707D'
orange = '#D55E00'
blue = '#0072B2'
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 9,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'axes.linewidth': 0.75,
        'axes.edgecolor': 'black',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': ink,
        'axes.labelcolor': ink,
        'xtick.color': 'black',
        'ytick.color': 'black',
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig, ax = plt.subplots(figsize = (3, 4), dpi = 150)
fig.subplots_adjust(left = 0.22, right = 0.97, bottom = 0.155, top = 0.935)
fig.text(0.035, 0.98, 'A', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    'Differential abundance of Input RNA',
    ha = 'left',
    va = 'center',
    fontsize = 11,
    fontweight = 'normal',
)
ax.set_axisbelow(True)
ax.grid(False)
y = -np.log10(df.adjustedPValue.clip(lower = 1e-300))
valid = df.adjustedPValue.notna() & np.isfinite(df.log2FoldChange)
bg = valid & ~sig
ax.scatter(
    df.loc[bg, 'log2FoldChange'],
    y[bg],
    s = 3.5,
    c = '#BCC4CD',
    alpha = 0.32,
    linewidths = 0,
    rasterized = True,
    zorder = 1,
)
for pos, color in [(False, blue), (True, orange)]:
    keep = sig & (df.log2FoldChange.gt(0) if pos else df.log2FoldChange.lt(0))
    ax.scatter(
        df.loc[keep, 'log2FoldChange'],
        y[keep],
        s = 25,
        c = color,
        edgecolors = 'black',
        linewidths = 0.5,
        zorder = 3,
    )
threshold = -np.log10(0.05)
ax.axhline(
    threshold, color = muted, linewidth = 0.85, linestyle = (0, (3, 2)), zorder = 0
)
threshold_text = ax.text(
    -13,
    threshold + 0.45,
    r'$P_{\mathrm{BH}}=0.05$',
    ha = 'center',
    va = 'bottom',
    fontsize = 7,
    color = muted,
    bbox = {'facecolor': 'white', 'edgecolor': 'none', 'pad': 1},
    zorder = 6,
)
ax.axvline(0, color = muted, linewidth = 0.7, zorder = 0)
genes = df.set_index('geneName')
gene_labels = []
for name, position in {'PHOX2A': (-11.0, 16.5), 'TOP2A': (6.0, 17.2)}.items():
    row = genes.loc[name]
    yy = -np.log10(row.adjustedPValue)
    ax.scatter(
        [row.log2FoldChange],
        [yy],
        s = 48,
        c = (orange if row.log2FoldChange > 0 else blue)
        if row.adjustedPValue < 0.05
        else '#BCC4CD',
        edgecolors = 'black',
        linewidths = 0.6,
        zorder = 4,
    )
    gene_labels.append(
        ax.annotate(
            name,
            (row.log2FoldChange, yy),
            xytext = position,
            fontsize = 10,
            fontstyle = 'italic',
            ha = 'center',
            va = 'center',
            bbox = {'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.5},
            arrowprops = {
                'arrowstyle': '-',
                'color': muted,
                'linewidth': 0.9,
                'shrinkA': 2,
                'shrinkB': 4,
            },
            zorder = 5,
        )
    )
ax.set(xlim = (-16.5, 10.5), ylim = (-0.8, 36.4))
ax.set_xticks([-15, -10, -5, 0, 5, 10])
ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel(r'Input log$_2$(ASD/Unaffected)', labelpad = 7)
ax.set_ylabel(r'$-\log_{10}(P_{\mathrm{BH}})$', labelpad = 7)
handles = [
    Line2D(
        [],
        [],
        linestyle = 'none',
        marker = 'o',
        markersize = 5,
        color = c,
        markeredgecolor = 'black',
        markeredgewidth = 0.5,
        label = label,
    )
    for c, label in [(orange, 'Higher in ASD (135)'), (blue, 'Lower in ASD (170)')]
]
ax.legend(
    handles = handles,
    loc = 'upper left',
    bbox_to_anchor = (-0.025, 1),
    fontsize = 8,
    frameon = False,
    handletextpad = 0.3,
    labelspacing = 0.4,
    borderpad = 0.35,
)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
points = ax.transData.transform(
    np.column_stack([df.loc[valid, 'log2FoldChange'], y[valid]])
)
radii = (
    np.where(sig[valid], np.sqrt(25) / 2 + 0.5, np.sqrt(3.5) / 2) * fig.dpi / 72
)
for item in [threshold_text, *gene_labels]:
    box = item.get_bbox_patch().get_window_extent(renderer)
    collisions = (
        (points[:, 0] + radii > box.x0)
        & (points[:, 0] - radii < box.x1)
        & (points[:, 1] + radii > box.y0)
        & (points[:, 1] - radii < box.y1)
    )
    assert not collisions.any(), (
        'Label background overlaps a gene marker',
        item.get_text(),
        int(collisions.sum()),
    )
outside = []
for t in fig.findobj(Text):
    if not t.get_visible() or not t.get_text().strip():
        continue
    b = t.get_window_extent(renderer)
    if b.x0 < 0 or b.y0 < 0 or b.x1 > fig.bbox.width or b.y1 > fig.bbox.height:
        outside.append(t.get_text())
assert not outside, outside
stem = 'Figure_2A_Input_volcano_titled_3x4'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (1800, 2400)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((750, 1000))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (216, 288))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 216, height = 288)
c.showPage()
c.save()
