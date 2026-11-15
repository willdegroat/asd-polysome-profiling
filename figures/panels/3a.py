import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('3A')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '3A'

P = PROJECT_ROOT
sources = ['transcriptGeneEvidence.tsv', 'transcriptGeneSensitivity.tsv']
hashes = {
    s: hashlib.sha256((P / 'results' / s).read_bytes()).hexdigest()
    for s in sources
}
d = pd.read_csv(P / 'results' / sources[0], sep = '\t')
d = d[d.analysis.eq('redistribution') & d.variant.eq('stable')].copy()
assert (
    len(d) == 11971
    and d.parentGeneId.is_unique
    and d.parentGeneAdjustedPValue.notna().all()
)
d['x'] = -np.log10(d.parentGeneAdjustedPValue)
d['y'] = -np.log10(d.geneAdjustedPValue)
d['category'] = np.select(
    [
        d.geneFdr & ~d.parentGeneFdr,
        d.geneFdr & d.parentGeneFdr,
        ~d.geneFdr & d.parentGeneFdr,
    ],
    ['Transcript-derived only', 'Both analyses', 'Gene-level only'],
    default = 'Neither',
)
order = [
    'Transcript-derived only',
    'Both analyses',
    'Gene-level only',
    'Neither',
]
counts = d.category.value_counts()
assert [int(counts.get(k, 0)) for k in order] == [48, 5, 9, 11909]
assert np.isfinite(d[['x', 'y']]).all().all()
ink = '#202938'
gray = '#63707D'
purple = '#6A51A3'
teal = '#D4B000'
colors = {
    'Transcript-derived only': purple,
    'Both analyses': teal,
    'Gene-level only': '#00B8D4',
    'Neither': '#C2C7CB',
}
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 9.5,
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
fig.text(0.035, 0.963, 'A', fontsize = 14, fontweight = 'bold', va = 'center')
title = fig.text(
    0.14,
    0.963,
    'Transcript analysis detects P/M changes in 48 additional genes',
    fontsize = 11,
    va = 'center',
)
ax = fig.add_axes([0.18, 0.19, 0.775, 0.715])
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = '#E7EBEF', lw = 0.55)
t = -np.log10(0.05)
ax.axvline(t, color = '#818B95', lw = 0.8, ls = (0, (3, 3)), zorder = 1)
ax.axhline(t, color = '#818B95', lw = 0.8, ls = (0, (3, 3)), zorder = 1)
for category in [
    'Neither',
    'Gene-level only',
    'Both analyses',
    'Transcript-derived only',
]:
    a = d[d.category.eq(category)]
    ax.scatter(
        a.x,
        a.y,
        s = 6 if category == 'Neither' else 26,
        c = colors[category],
        marker = 'o',
        alpha = 0.25 if category == 'Neither' else 1,
        edgecolors = 'none' if category == 'Neither' else 'black',
        linewidths = 0 if category == 'Neither' else 0.55,
        zorder = 2 if category == 'Neither' else 3,
        rasterized = category == 'Neither',
    )
ax.set_xlim(-0.06, 4.05)
ax.set_ylim(-0.12, 8.15)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_yticks([0, 2, 4, 6, 8])
ax.tick_params(labelsize = 8.5)
ax.set_xlabel(
    'Gene-level P/M evidence\n' + r'$-\log_{10}(P_{\mathrm{BH}})$', labelpad = 7
)
ax.set_ylabel(
    'Transcript-derived P/M evidence\n' + r'$-\log_{10}(P_{\mathrm{BH}})$',
    labelpad = 7,
)
ax.text(
    t + 0.055,
    7.90,
    r'$P_{\mathrm{BH}}=0.05$',
    rotation = 90,
    ha = 'left',
    va = 'top',
    fontsize = 7.5,
    color = gray,
)
ax.text(
    3.96,
    t + 0.08,
    r'$P_{\mathrm{BH}}=0.05$',
    ha = 'right',
    va = 'bottom',
    fontsize = 7.5,
    color = gray,
)
labels = {
    'DLK1': (0.62, 3.30),
    'COPS7A': (0.06, 7.85),
    'NEUROG1': (2.35, 6.20),
    'TOP2A': (2.10, 0.80),
    'PHOX2A': (3.00, 0.55),
}
for gene, position in labels.items():
    r = d[d.geneName.eq(gene)].iloc[0]
    assert (r.x - t) * (position[0] - t) > 0 and (r.y - t) * (
        position[1] - t
    ) > 0
    ax.annotate(
        gene,
        xy = (r.x, r.y),
        xytext = position,
        fontsize = 8.5,
        fontstyle = 'italic',
        ha = 'left',
        va = 'center',
        arrowprops = {
            'arrowstyle': '-',
            'color': '#89939B',
            'lw': 0.75,
            'shrinkA': 3,
            'shrinkB': 4,
        },
        zorder = 4,
    )
category_labels = {
    'Transcript-derived only': 'Transcript only',
    'Both analyses': 'Both',
    'Gene-level only': 'Gene only',
    'Neither': 'Neither',
}
handles = [
    Line2D(
        [],
        [],
        marker = 'o',
        color = colors[k],
        linestyle = 'none',
        markeredgecolor = 'none' if k == 'Neither' else 'black',
        markeredgewidth = 0.5,
        markersize = 5,
        label = f'{category_labels[k]} ({counts[k]:,})',
    )
    for k in [
        'Transcript-derived only',
        'Gene-level only',
        'Both analyses',
        'Neither',
    ]
]
legend = fig.legend(
    handles = handles,
    loc = 'lower center',
    bbox_to_anchor = (0.5675, 0.010),
    ncol = 2,
    frameon = False,
    fontsize = 7.8,
    handlelength = 0.8,
    columnspacing = 0.9,
    handletextpad = 0.35,
    labelspacing = 0.35,
    borderpad = 0.1,
)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
while title.get_window_extent(renderer).x1 > fig.bbox.width * 0.975:
    title.set_fontsize(title.get_fontsize() - 0.2)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
outside = []
for text in fig.findobj(Text):
    if text.get_visible() and text.get_text().strip():
        box = text.get_window_extent(renderer)
        if (
            box.x0 < 0
            or box.y0 < 0
            or box.x1 > fig.bbox.width
            or box.y1 > fig.bbox.height
        ):
            outside.append(text.get_text())
assert not outside, outside
assert not legend.get_window_extent(renderer).overlaps(
    ax.xaxis.label.get_window_extent(renderer)
)
assert d.x.max() < ax.get_xlim()[1] and d.y.max() < ax.get_ylim()[1]
stem = 'Figure_3A_transcript_vs_gene_PM_evidence_clean_5x6'
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
