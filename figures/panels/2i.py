import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.text import Text
from matplotlib.transforms import blended_transform_factory
from PIL import Image, ImageOps
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2I')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2I'

"""SynGO + Fu dot-only Figure 3E; 5 x 4 in, 600 DPI."""
SOURCE = PROJECT_ROOT / 'results/rankedEnrichment.tsv'
C = {
    'ink': '#202938',
    'muted': '#63707D',
    'grid': '#E7EBEF',
    'monosome': '#526779',
    'polysome': '#B5367A',
}
source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
table = pd.read_csv(SOURCE, sep = '\t')
data = table[
    table.analysis.eq('redistribution')
    & table.collection.isin(['synGO', 'fuAsd'])
    & table.adjustedPValue.lt(0.05)
].copy()
assert len(data) == 9 and data.collection.value_counts().to_dict() == {
    'synGO': 7,
    'fuAsd': 2,
}
data = data.sort_values(
    'normalizedEnrichmentScore', ascending = False
).reset_index(drop = True)
names = {
    (
        'Anchored Component Of Synaptic Vesicle Membrane (G'
        'O:0098993) CC'
    ): (
        'Vesicle membrane (anchored)'
    ),
    (
        'Integral Component Of Synaptic Vesicle Membrane (G'
        'O:0030285) CC'
    ): (
        'Vesicle membrane (integral)'
    ),
    'Synaptic Vesicle Membrane (GO:0030672) CC': 'Synaptic vesicle membrane',
    'Postsynaptic Ribosome CC': 'Postsynaptic ribosome',
    'Presynaptic Ribosome CC': 'Presynaptic ribosome',
    'Translation At Presynapse (GO:0140236) BP': 'Translation at presynapse',
    'Translation At Postsynapse (GO:0140242) BP': 'Translation at postsynapse',
    'Fu ASD FDR <= 0.001': 'ASD72',
    'Fu ASD FDR <= 0.05': 'ASD185',
}
data['displayLabel'] = [
    names[r.term]
    + (r' [Fu $\it{et\ al.}$]' if r.collection == 'fuAsd' else ' [SynGO]')
    for r in data.itertuples()
]


def size(q):
    return 16 + 32 * (-np.log10(np.asarray(q)) - (-np.log10(0.05)))


data['dotAreaPointsSquared'] = size(data.adjustedPValue)
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'text.color': C['ink'],
        'axes.labelcolor': C['ink'],
        'axes.edgecolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.linewidth': 0.7,
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'savefig.facecolor': 'white',
        'figure.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (5, 4), dpi = 150)
fig.text(0.035, 0.98, 'I', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    'ASD-associated gene enrichment supports\na shift toward monosomes in ASD',
    ha = 'left',
    va = 'center',
    multialignment = 'left',
    fontsize = 9.25,
    fontweight = 'normal',
    linespacing = 1.0,
)
ax = fig.add_axes([0.48, 0.29, 0.49, 0.57])
ax.set(xlim = (-3, 3), ylim = (8.5, -0.5))
ax.set_yticks(np.arange(9))
ax.tick_params(axis = 'y', length = 0, labelleft = False)
ax.set_xticks(
    [-3, -2, -1, 0, 1, 2, 3], ['−3', '−2', '−1', '0', '1', '2', '3'], fontsize = 9
)
ax.tick_params(axis = 'x', length = 3, width = 0.6, pad = 3)
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = C['grid'], linewidth = 0.6)
for key in ['left', 'top', 'right']:
    ax.spines[key].set_visible(False)
ax.set_xlabel('NES', fontsize = 9, labelpad = 5)
labels = []
transform = blended_transform_factory(fig.transFigure, ax.transData)
for i, row in data.iterrows():
    labels.append(
        ax.text(
            0.04,
            i,
            row.displayLabel,
            transform = transform,
            ha = 'left',
            va = 'center',
            fontsize = 8.5,
            clip_on = False,
        )
    )
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
right = max(t.get_window_extent(renderer).x1 for t in labels)
left = (right + 7 * fig.dpi / 72) / fig.bbox.width
assert left < 0.56, left
ax.set_position([left, 0.29, 0.97 - left, 0.57])
cmap = LinearSegmentedColormap.from_list(
    'Monopoly_NES', [C['monosome'], '#F7F7F7', C['polysome']]
)
norm = TwoSlopeNorm(vmin = -3, vcenter = 0, vmax = 3)
points = ax.scatter(
    data.normalizedEnrichmentScore,
    np.arange(9),
    s = data.dotAreaPointsSquared,
    c = data.normalizedEnrichmentScore,
    cmap = cmap,
    norm = norm,
    edgecolors = 'black',
    linewidths = 0.45,
    zorder = 3,
)
legend_q = [0.001, 0.005, 0.02]
handles = [
    ax.scatter(
        [],
        [],
        s = float(size(q)),
        facecolors = C['muted'],
        edgecolors = 'black',
        linewidths = 0.45,
        label = f'{q:g}',
    )
    for q in legend_q
]
legend = fig.legend(
    handles = handles,
    title = 'FDR ' + r'$q$',
    loc = 'lower center',
    bbox_to_anchor = ((left + 0.97) / 2, 0.015),
    ncol = 3,
    fontsize = 8,
    title_fontsize = 8,
    frameon = False,
    handlelength = 1.2,
    handletextpad = 0.45,
    columnspacing = 1.1,
    borderpad = 0.25,
    labelspacing = 0.45,
)
interpretation = []
for yy, heading, finding, color in [
    (0.225, 'ASD-associated genes:', 'Lower P/M in ASD', '#0072B2'),
    (0.115, 'Synaptic vesicle genes:', 'Higher P/M in ASD', '#D55E00'),
]:
    interpretation.append(
        fig.text(
            0.04, yy, heading, ha = 'left', va = 'top', fontsize = 8.5, color = C['ink']
        )
    )
    interpretation.append(
        fig.text(
            0.04,
            yy - 0.043,
            finding,
            ha = 'left',
            va = 'top',
            fontsize = 8.5,
            color = color,
            fontweight = 'bold',
        )
    )
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
boxes = [t.get_window_extent(renderer) for t in labels]
assert all(
    not a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1 :]
)
assert not legend.get_window_extent(renderer).overlaps(
    ax.xaxis.label.get_window_extent(renderer)
), 'Legend/axis label overlap.'
for t in interpretation:
    assert not t.get_window_extent(renderer).overlaps(
        legend.get_window_extent(renderer)
    ), 'Interpretation/legend overlap.'
outside = []
for text in fig.findobj(Text):
    if not text.get_visible() or not text.get_text().strip():
        continue
    b = text.get_window_extent(renderer)
    if (
        b.x0 < -0.3
        or b.y0 < -0.3
        or b.x1 > fig.bbox.width + 0.3
        or b.y1 > fig.bbox.height + 0.3
    ):
        outside.append(text.get_text())
assert not outside, outside
assert len(ax.lines) == 0
stem = 'Figure_2I_SynGO_Fu_dotplot_titled_5x4'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 288))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 288)
c.showPage()
c.save()
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 2400) and np.allclose(
        im.info['dpi'], [600, 600], atol = 0.1
    )
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    ImageOps.contain(im.convert('RGB'), (1200, 960)).save(REV / 'preview.png')
