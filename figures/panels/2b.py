import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import blended_transform_factory
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2B')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2B'

source = PROJECT_ROOT / 'results/overrepresentationAnalysis.tsv'
digest = hashlib.sha256(source.read_bytes()).hexdigest()
table = pd.read_csv(source, sep = '\t')
eligible = table[
    table.analysis.eq('input')
    & table.collection.isin(
        ['goBiologicalProcess', 'goCellularComponent', 'goMolecularFunction']
    )
    & table.adjustedPValue.lt(0.05)
]
data = (
    eligible.sort_values(
        ['adjustedPValue', 'pValue', 'collection', 'term', 'direction'],
        kind = 'stable',
    )
    .head(6)
    .copy()
    .reset_index(drop = True)
)
assert len(data) == 6 and data.direction.eq('negative').all()
names = {
    'Nervous System Development (GO:0007399)': 'Nervous system\ndevelopment [GO BP]',
    (
        'Regulation of Osteoclast Differentiation (GO:00456'
        '70)'
    ): (
        'Regulation of osteoclast\ndifferentiation [GO BP]'
    ),
    (
        'DNA-binding Transcription Activator Activity, RNA '
        'Polymerase II-specific (GO:0001228)'
    ): (
        'DNA-binding transcription\nactivator activity,\nRNA '
        'Pol II-specific [GO MF]'
    ),
    (
        'Positive Regulation of Cell Population Proliferati'
        'on (GO:0008284)'
    ): (
        'Positive regulation of\ncell population proliferati'
        'on\n[GO BP]'
    ),
    'Neuron Differentiation (GO:0030182)': 'Neuron differentiation\n[GO BP]',
    'Generation of Neurons (GO:0048699)': 'Generation of neurons\n[GO BP]',
}
data['displayLabel'] = data.term.map(names)
assert data.displayLabel.notna().all()
neural_terms = {
    'Nervous System Development (GO:0007399)',
    'Neuron Differentiation (GO:0030182)',
    'Generation of Neurons (GO:0048699)',
}
data['neuralTermLabelBold'] = data.term.isin(neural_terms)
data['dotAreaPointsSquared'] = 20 + 6 * data.overlapGenesCount
ink = '#202938'
muted = '#63707D'
blue = '#0072B2'
orange = '#D55E00'
grid = '#E7EBEF'
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 9,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'text.color': ink,
        'axes.labelcolor': ink,
        'axes.edgecolor': 'black',
        'axes.linewidth': 0.75,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.color': 'black',
        'ytick.color': 'black',
        'xtick.labelsize': 8.5,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (4, 4), dpi = 150)
fig.text(0.035, 0.98, 'B', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    'Genes with lower Input abundance\nare enriched for neuronal development',
    ha = 'left',
    va = 'center',
    multialignment = 'left',
    fontsize = 9.25,
    fontweight = 'normal',
    linespacing = 1.0,
)
ax = fig.add_axes([0.54, 0.23, 0.43, 0.69])
ax.set(xlim = (0, 21), ylim = (5.5, -0.5))
ax.set_yticks(range(6))
ax.tick_params(axis = 'y', length = 0, labelleft = False)
ax.set_xticks([0, 5, 10, 15, 20])
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = grid, linewidth = 0.65)
ax.spines['left'].set_visible(False)
ax.set_xlabel('Fold enrichment', labelpad = 5)
transform = blended_transform_factory(fig.transFigure, ax.transData)
labels = [
    ax.text(
        0.035,
        i,
        row.displayLabel,
        transform = transform,
        ha = 'left',
        va = 'center',
        fontsize = 8.5,
        fontweight = 'bold' if row.neuralTermLabelBold else 'normal',
        linespacing = 1.15,
        clip_on = False,
    )
    for i, row in data.iterrows()
]
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
left = (
    max(t.get_window_extent(renderer).x1 for t in labels) + 7 * fig.dpi / 72
) / fig.bbox.width
assert left < 0.64, left
ax.set_position([left, 0.23, 0.97 - left, 0.69])
ax.scatter(
    data.enrichmentRatio,
    np.arange(6),
    s = data.dotAreaPointsSquared,
    c = [blue if d == 'negative' else orange for d in data.direction],
    edgecolors = 'black',
    linewidths = 0.55,
    zorder = 3,
)
direction_handles = [
    Line2D(
        [],
        [],
        marker = 'o',
        linestyle = 'none',
        color = c,
        markeredgecolor = 'black',
        markeredgewidth = 0.55,
        markersize = 5,
        label = label,
    )
    for c, label in [(orange, 'Higher in ASD'), (blue, 'Lower in ASD')]
]
direction_legend = fig.legend(
    handles = direction_handles,
    title = 'Input gene list:',
    loc = 'lower left',
    bbox_to_anchor = (0.035, 0.015),
    alignment = 'left',
    frameon = False,
    fontsize = 8,
    title_fontsize = 8,
    handletextpad = 0.5,
    handlelength = 1,
    labelspacing = 0.4,
    borderpad = 0.2,
)
counts = [5, 10, 20]
handles = [
    ax.scatter(
        [],
        [],
        s = 20 + 6 * n,
        c = muted,
        edgecolors = 'black',
        linewidths = 0.55,
        label = str(n),
    )
    for n in counts
]
legend = fig.legend(
    handles = handles,
    title = 'Overlapping genes',
    loc = 'lower center',
    bbox_to_anchor = ((left + 0.97) / 2, 0.025),
    ncol = 3,
    frameon = False,
    fontsize = 8,
    title_fontsize = 8,
    handletextpad = 0.35,
    handlelength = 1.1,
    columnspacing = 0.8,
    labelspacing = 0.55,
    borderpad = 0.2,
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
boxes = [t.get_window_extent(renderer) for t in labels]
assert all(
    not a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1 :]
)
assert not legend.get_window_extent(renderer).overlaps(
    ax.xaxis.label.get_window_extent(renderer)
)
assert not direction_legend.get_window_extent(renderer).overlaps(
    legend.get_window_extent(renderer)
)
stem = 'Figure_2B_Input_GO_ORA_titled_4x4'
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
