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

OUT = panel_output('2D')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2D'

ROOT = HERE.parent
source = PROJECT_ROOT / 'results/polysomeMonosome.tsv'
before = hashlib.sha256(source.read_bytes()).hexdigest()
df = pd.read_csv(source, sep = '\t')
sig = df.adjustedPValue.lt(0.05)
assert len(df) == 20470 and sig.sum() == 22
assert (sig & df.log2FoldChange.gt(0)).sum() == 12
assert (sig & df.log2FoldChange.lt(0)).sum() == 10
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
fig, ax = plt.subplots(figsize = (4, 4), dpi = 150)
fig.subplots_adjust(left = 0.16, right = 0.97, bottom = 0.155, top = 0.935)
fig.text(0.035, 0.98, 'D', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    'Differential P/M ratios',
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
ax.text(
    0.98,
    threshold + 0.08,
    r'$P_{\mathrm{BH}}=0.05$',
    transform = ax.get_yaxis_transform(),
    ha = 'right',
    va = 'bottom',
    fontsize = 8,
    color = muted,
    bbox = {'facecolor': 'white', 'edgecolor': 'none', 'pad': 2},
    zorder = 6,
)
ax.axvline(0, color = muted, linewidth = 0.7, zorder = 0)
genes = df.set_index('geneName')
for name, position in {'PHOX2A': (3.4, 4.02), 'TOP2A': (-5.1, 3.12)}.items():
    row = genes.loc[name]
    yy = -np.log10(row.adjustedPValue)
    ax.scatter(
        [row.log2FoldChange],
        [yy],
        s = 48,
        c = orange if row.log2FoldChange > 0 else blue,
        edgecolors = 'black',
        linewidths = 0.6,
        zorder = 4,
    )
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
ax.set(xlim = (-11.5, 17), ylim = (-0.06, 4.3))
ax.set_xticks([-10, -5, 0, 5, 10, 15])
ax.set_yticks([0, 1, 2, 3, 4])
ax.set_xlabel(r'Change in P/M (log$_2$ ASD/Unaffected)', labelpad = 7)
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
    for c, label in [(orange, 'Higher P/M (12)'), (blue, 'Lower P/M (10)')]
]
ax.legend(
    handles = handles,
    loc = 'upper right',
    bbox_to_anchor = (1, 1),
    fontsize = 8,
    frameon = True,
    facecolor = 'white',
    edgecolor = 'none',
    framealpha = 1,
    handletextpad = 0.3,
    labelspacing = 0.4,
    borderpad = 0.35,
)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
outside = []
for t in fig.findobj(Text):
    if not t.get_visible() or not t.get_text().strip():
        continue
    b = t.get_window_extent(renderer)
    if b.x0 < 0 or b.y0 < 0 or b.x1 > fig.bbox.width or b.y1 > fig.bbox.height:
        outside.append(t.get_text())
assert not outside, outside
stem = 'Figure_2D_PM_volcano_titled_4x4'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (2400, 2400)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1000, 1000))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (288, 288))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 288, height = 288)
c.showPage()
c.save()
