import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2H')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2H'

PROJECT = PROJECT_ROOT
source = PROJECT / 'results/fuAsdLeaveOneOutEnrichment.tsv'
primary_source = PROJECT / 'results/fuAsdRankedEnrichment.tsv'
metadata_source = PROJECT / 'data/metadata.tsv'
sources = [
    {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    for p in [source, primary_source, metadata_source]
]
loo = pd.read_csv(source, sep = '\t', dtype = {'omittedParticipantId': str})
primary = pd.read_csv(primary_source, sep = '\t')
meta = pd.read_csv(metadata_source, sep = '\t', dtype = str)
meta['displayLabel'] = [
    f'U{i + 1}' if r.condition == 'Unaffected' else f'ASD{i - 2}'
    for i, r in enumerate(meta.itertuples())
]
assert meta.displayLabel.tolist() == ['U1', 'U2', 'U3', 'ASD1', 'ASD2', 'ASD3']
mapping = meta[['sampleId', 'condition', 'displayLabel']].copy()
tiers = {
    'Fu ASD FDR \u2264 0.001': ('ASD72', 'o', -0.065),
    'Fu ASD FDR \u2264 0.05': ('ASD185', 'D', 0.065),
}
loo['displayLabel'] = loo.omittedParticipantId.map(
    mapping.set_index('sampleId').displayLabel
)
loo['geneSet'] = loo.tier.map({k: v[0] for k, v in tiers.items()})
loo['significant'] = loo.tierAdjustedPValue.lt(0.05)
export = loo[
    [
        'displayLabel',
        'omittedParticipantId',
        'omittedCondition',
        'geneSet',
        'normalizedEnrichmentScore',
        'tierAdjustedPValue',
        'significant',
    ]
]
assert (
    len(export) == 12
    and export.groupby('geneSet').significant.sum().eq(5).all()
)
assert loo.normalizedEnrichmentScore.abs().max() <= 3
future = pd.read_csv(PROJECT / 'results/rankedEnrichment.tsv', sep = '\t')
future = future[
    future.analysis.eq('redistribution')
    & future.collection.isin(['synGO', 'fuAsd'])
    & future.adjustedPValue.lt(0.05)
]
assert len(future) == 9 and future.normalizedEnrichmentScore.abs().max() <= 3
ink = '#202938'
muted = '#63707D'
palette = {'negative': '#526779', 'zero': '#F7F7F7', 'positive': '#B5367A'}
cmap = LinearSegmentedColormap.from_list(
    'NES_slate_white_magenta', list(palette.values())
)
norm = TwoSlopeNorm(vmin = -3, vcenter = 0, vmax = 3)
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 8.5,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'text.color': ink,
        'axes.labelcolor': ink,
        'axes.edgecolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.linewidth': 0.75,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (5, 2), dpi = 150)
fig.text(0.02, 0.98, 'H', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.945,
    'ASD-associated enrichment is sensitive to ASD1 omission',
    ha = 'left',
    va = 'center',
    fontsize = 11,
    fontweight = 'normal',
)
ax = fig.add_axes([0.14, 0.27, 0.60, 0.60])
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = '#E7EBEF', linewidth = 0.55)
baselines = []
for tier, (name, marker, offset) in tiers.items():
    frame = (
        loo[loo.tier.eq(tier)]
        .set_index('omittedParticipantId')
        .loc[meta.sampleId]
    )
    x = np.arange(6) + offset
    y = frame.normalizedEnrichmentScore.to_numpy()
    # Color interpolates continuously with the NES along each segment.
    xx = np.concatenate(
        [np.linspace(x[i], x[i + 1], 65, endpoint = False) for i in range(5)]
        + [x[-1:]]
    )
    yy = np.interp(xx, x, y)
    points = np.column_stack([xx, yy])
    segments = np.stack([points[:-1], points[1:]], axis = 1)
    line = LineCollection(
        segments, cmap = cmap, norm = norm, linewidths = 1.1, zorder = 2
    )
    line.set_array((yy[:-1] + yy[1:]) / 2)
    ax.add_collection(line)
    for xpos, value, significant in zip(
        x, y, frame.tierAdjustedPValue.lt(0.05)
    ):
        ax.plot(
            xpos,
            value,
            marker = marker,
            markersize = 4.7,
            linestyle = 'none',
            markeredgewidth = 0.6,
            markeredgecolor = 'black',
            markerfacecolor = cmap(norm(value)) if significant else 'white',
            zorder = 4,
        )
    value = float(
        primary[
            (primary.context == 'redistribution') & primary.tier.eq(tier)
        ].normalizedEnrichmentScore.iloc[0]
    )
    baselines.append({'geneSet': name, 'fullDataNES': value})
ax.set(xlim = (-0.4, 5.4), ylim = (-2.85, 1.35))
ax.set_xticks(range(6), meta.displayLabel, fontsize = 8.5)
ax.set_yticks([-2, -1, 0, 1], ['\u22122', '\u22121', '0', '1'], fontsize = 8.5)
ax.tick_params(length = 2.4, width = 0.6, pad = 2.5)
ax.set_ylabel('NES', fontsize = 9, labelpad = 5)
ax.set_xlabel('Sample omitted', fontsize = 8.5, labelpad = 4)
handles = [
    Line2D(
        [],
        [],
        marker = 'o',
        color = ink,
        linestyle = 'none',
        markersize = 4.7,
        label = 'ASD72',
    ),
    Line2D(
        [],
        [],
        marker = 'D',
        color = ink,
        linestyle = 'none',
        markersize = 4.7,
        label = 'ASD185',
    ),
    tuple(
        Line2D(
            [],
            [],
            marker = m,
            markerfacecolor = 'white',
            markeredgecolor = 'black',
            linestyle = 'none',
            markersize = 4.7,
        )
        for m in ['o', 'D']
    ),
]
fig.legend(
    handles = handles,
    labels = ['ASD72', 'ASD185', r'$P_{\mathrm{BH}}\geq0.05$'],
    handler_map = {tuple: HandlerTuple(ndivide = None, pad = 0.4)},
    loc = 'upper left',
    bbox_to_anchor = (0.75, 0.89),
    fontsize = 7.8,
    frameon = False,
    handlelength = 1.65,
    handletextpad = 0.5,
    labelspacing = 0.7,
    borderpad = 0.25,
    borderaxespad = 0,
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
stem = 'Figure_2H_sample_omission_titled_5x2'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 1200)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1500, 600))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 144))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 144)
c.showPage()
c.save()
