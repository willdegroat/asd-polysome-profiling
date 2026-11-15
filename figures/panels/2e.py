import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2E')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2E'

profile_path = HERE / 'PHOX2A_source_samples.csv'
pm_path = PROJECT_ROOT / 'results/polysomeMonosome.tsv'
sources = [
    {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    for p in [profile_path, pm_path]
]
data = pd.read_csv(profile_path)
data = data[data.geneName.eq('PHOX2A')].copy()
assert (
    len(data) == 18
    and data.groupby(['condition', 'fraction']).size().eq(3).all()
)
assert np.allclose(data.log2TPMPlus1, np.log2(data.TPM + 1))
means = data.groupby(['condition', 'fraction'], as_index = False).agg(
    meanLog2TPMPlus1 = ('log2TPMPlus1', 'mean'), n = ('sampleId', 'size')
)
row = pd.read_csv(pm_path, sep = '\t').set_index('geneName').loc['PHOX2A']
fold = float(np.exp2(row.log2FoldChange))
ink = '#202938'
muted = '#63707D'
orange = '#D55E00'
colors = {'Unaffected': '#238B45', 'ASD': '#C43C39'}
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
fig.subplots_adjust(left = 0.22, right = 0.97, bottom = 0.19, top = 0.795)
fig.text(0.035, 0.98, 'E', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    r'$\it{PHOX2A}$ abundance across RNA fractions',
    ha = 'left',
    va = 'center',
    fontsize = 9.5,
    fontweight = 'normal',
)
center = (0.22 + 0.97) / 2
fig.text(
    center,
    0.93,
    'PHOX2A',
    ha = 'center',
    va = 'top',
    fontsize = 11,
    fontstyle = 'italic',
    fontweight = 'bold',
)
fig.text(
    center,
    0.888,
    'Fitted P/M ratio',
    ha = 'center',
    va = 'top',
    fontsize = 8.5,
    color = muted,
)
fig.text(
    center,
    0.853,
    f'{fold:.1f}-fold higher in ASD',
    ha = 'center',
    va = 'top',
    fontsize = 9,
    color = orange,
    fontweight = 'bold',
)
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = '#E7EBEF', linewidth = 0.65)
fraction_markers = {'input': 'o', 'monosome': '^', 'polysome': 's'}
plotted_points = []
for condition in ['Unaffected', 'ASD']:
    values = (
        means[means.condition.eq(condition)]
        .set_index('fraction')
        .loc[['input', 'monosome', 'polysome'], 'meanLog2TPMPlus1']
    )
    ax.plot(
        range(3),
        values,
        color = colors[condition],
        linewidth = 1.8,
        linestyle = '-',
        zorder = 3,
    )
    for xpos, (fraction, value) in enumerate(values.items()):
        ax.plot(
            xpos,
            value,
            color = colors[condition],
            marker = fraction_markers[fraction],
            markersize = 5.6,
            markeredgewidth = 0.6,
            markeredgecolor = 'black',
            linestyle = 'none',
            zorder = 3,
        )
        plotted_points.append(
            {
                'condition': condition,
                'fraction': fraction,
                'x': xpos,
                'y': float(value),
                'marker': fraction_markers[fraction],
            }
        )
assert len(plotted_points) == 6
ax.set(ylim = (-0.12, 7), xlim = (-0.28, 2.28))
ax.set_xticks(range(3), ['Input', 'Monosome', 'Polysome'])
ax.set_yticks([0, 2, 4, 6])
ax.set_ylabel(r'$\log_{2}(\mathrm{TPM}+1)$', labelpad = 7)
handles = [
    tuple(
        Line2D(
            [],
            [],
            color = colors[c],
            marker = m,
            markeredgecolor = 'black',
            markeredgewidth = 0.6,
            linestyle = 'none',
            markersize = 5.6,
        )
        for m in ['o', '^', 's']
    )
    for c in ['Unaffected', 'ASD']
]
fig.legend(
    handles = handles,
    labels = ['Unaffected', 'ASD'],
    handler_map = {tuple: HandlerTuple(ndivide = None, pad = 0.25)},
    loc = 'lower center',
    bbox_to_anchor = (0.56, 0.025),
    ncol = 2,
    frameon = False,
    fontsize = 8.5,
    handlelength = 3.0,
    columnspacing = 1.2,
    handletextpad = 0.5,
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
labels = [t.get_window_extent(renderer) for t in ax.get_xticklabels()]
assert all(a.x1 < b.x0 for a, b in zip(labels, labels[1:])), (
    'Fraction labels overlap'
)
stem = 'Figure_2E_PHOX2A_fraction_shapes_titled_3x4'
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
