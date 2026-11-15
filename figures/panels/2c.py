import hashlib
import json

import matplotlib
import numpy as np
import pandas as pd
from _common import ASSETS, PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2C')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2C'

PROJECT = PROJECT_ROOT
paths = [
    PROJECT / 'results' / f
    for f in ['inputDifferentialExpression.tsv', 'polysomeMonosome.tsv']
]
hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
inp = pd.read_csv(paths[0], sep = '\t')
pm = pd.read_csv(paths[1], sep = '\t')
annotations = json.loads((ASSETS / 'gene_biotypes.json').read_text())['records']
significant = pm[pm.adjustedPValue.lt(0.05)].copy()
significant['geneType'] = significant.geneId.map(
    {k: v['geneType'] for k, v in annotations.items()}
)
significant['annotationGeneName'] = significant.geneId.map(
    {k: v['annotationGeneName'] for k, v in annotations.items()}
)
assert significant.geneName.eq(significant.annotationGeneName).all()
eligible = significant[significant.geneType.eq('protein_coding')].copy()
assert len(eligible) == 17
eligible['absolutePMLog2FoldChange'] = eligible.log2FoldChange.abs()
eligible = eligible.sort_values(
    ['absolutePMLog2FoldChange', 'geneId'],
    ascending = [False, True],
    kind = 'stable',
).reset_index(drop = True)
eligible['effectSizeRank'] = np.arange(1, len(eligible) + 1)
top = eligible.head(10).copy()
genes = top.geneName.tolist()
assert genes == [
    'H2BC17',
    'DHRSX',
    'DNAI3',
    'TOP2A',
    'NPIPB3',
    'PHOX2A',
    'NEUROG1',
    'SLC31A1',
    'GCN1',
    'ATP1B2',
]
columns = [
    'geneId',
    'geneName',
    'log2FoldChange',
    'log2FoldChangeSE',
    'pValue',
    'adjustedPValue',
]
selected = top[
    columns + ['geneType', 'absolutePMLog2FoldChange', 'effectSizeRank']
].merge(
    inp[columns],
    on = 'geneId',
    how = 'left',
    suffixes = ('PM', 'Input'),
    validate = 'one_to_one',
)
selected = selected.sort_values('effectSizeRank').reset_index(drop = True)
assert (
    len(selected) == 10 and selected.geneNamePM.eq(selected.geneNameInput).all()
)
assert (
    np.isfinite(
        selected[
            [
                'log2FoldChangeInput',
                'log2FoldChangePM',
                'log2FoldChangeSEInput',
                'log2FoldChangeSEPM',
                'adjustedPValueInput',
                'adjustedPValuePM',
            ]
        ]
    )
    .all()
    .all()
)
selected['InputSignificant'] = selected.adjustedPValueInput.lt(0.05)
selected['PMSignificant'] = selected.adjustedPValuePM.lt(0.05)
assert selected.PMSignificant.all()
assert set(selected.loc[selected.InputSignificant, 'geneNamePM']) == {
    'PHOX2A',
    'NEUROG1',
}
critical = 1.959963984540054
for contrast in ['Input', 'PM']:
    selected[f'ci95Lower{contrast}'] = (
        selected[f'log2FoldChange{contrast}']
        - critical * selected[f'log2FoldChangeSE{contrast}']
    )
    selected[f'ci95Upper{contrast}'] = (
        selected[f'log2FoldChange{contrast}']
        + critical * selected[f'log2FoldChangeSE{contrast}']
    )
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
fig, ax = plt.subplots(figsize = (3, 4), dpi = 150)
fig.subplots_adjust(left = 0.30, right = 0.97, bottom = 0.205, top = 0.935)
fig.text(0.035, 0.98, 'C', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.963,
    'P/M analysis reveals changes\nbeyond Input abundance',
    ha = 'left',
    va = 'center',
    multialignment = 'left',
    fontsize = 9.25,
    fontweight = 'normal',
    linespacing = 1.0,
)
yy = np.arange(len(selected))
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = grid, linewidth = 0.55)
ax.axvline(0, color = muted, linewidth = 0.7, zorder = 0)
for i, row in selected.iterrows():
    for contrast, marker, offset in [('Input', 'o', -0.15), ('PM', 'D', 0.15)]:
        estimate = row[f'log2FoldChange{contrast}']
        se = row[f'log2FoldChangeSE{contrast}']
        color = (
            muted if contrast == 'Input' else orange if estimate > 0 else blue
        )
        facecolor = (
            'white'
            if contrast == 'Input' and not row.InputSignificant
            else color
        )
        ax.errorbar(
            estimate,
            i + offset,
            xerr = critical * se,
            fmt = marker,
            markersize = 4,
            color = color,
            ecolor = color,
            elinewidth = 0.9,
            capthick = 0.7,
            capsize = 1.6,
            markerfacecolor = facecolor,
            markeredgecolor = 'black',
            markeredgewidth = 0.5,
            zorder = 3,
        )
ax.set_yticks(yy, [r'$\it{' + g + r'}$' for g in genes], fontsize = 8.2)
ax.tick_params(axis = 'y', length = 0, pad = 5)
ax.spines['left'].set_visible(False)
ax.set(ylim = (9.55, -0.55), xlim = (-8, 8))
ax.set_xticks([-8, -4, 0, 4, 8])
ax.set_xlabel(r'log$_2$(ASD/Unaffected)', labelpad = 6)
handles = [
    Line2D(
        [],
        [],
        marker = 'o',
        color = muted,
        linestyle = 'none',
        markersize = 4,
        markeredgecolor = 'black',
        markeredgewidth = 0.5,
        label = r'Input: $P_{\mathrm{BH}}<0.05$',
    ),
    Line2D(
        [],
        [],
        marker = 'o',
        color = muted,
        markerfacecolor = 'white',
        linestyle = 'none',
        markersize = 4,
        markeredgecolor = 'black',
        markeredgewidth = 0.5,
        label = r'Input: $P_{\mathrm{BH}}\geq0.05$',
    ),
    Line2D(
        [],
        [],
        marker = 'D',
        color = orange,
        linestyle = 'none',
        markersize = 4,
        markeredgecolor = 'black',
        markeredgewidth = 0.5,
        label = 'Higher P/M',
    ),
    Line2D(
        [],
        [],
        marker = 'D',
        color = blue,
        linestyle = 'none',
        markersize = 4,
        markeredgecolor = 'black',
        markeredgewidth = 0.5,
        label = 'Lower P/M',
    ),
]
legend = fig.legend(
    handles = handles,
    loc = 'lower center',
    bbox_to_anchor = (0.59, 0.020),
    alignment = 'left',
    ncol = 2,
    fontsize = 7.2,
    frameon = False,
    handlelength = 0.8,
    handletextpad = 0.3,
    columnspacing = 0.75,
    labelspacing = 0.45,
    borderpad = 0.2,
)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
outside = []
excluded = set()
for axis in [ax.xaxis, ax.yaxis]:
    lo, hi = sorted(axis.get_view_interval())
    for tick in axis.get_major_ticks():
        if not lo <= tick.get_loc() <= hi:
            excluded.update([id(tick.label1), id(tick.label2)])
for t in fig.findobj(Text):
    if not t.get_visible() or not t.get_text().strip() or id(t) in excluded:
        continue
    b = t.get_window_extent(renderer)
    if b.x0 < 0 or b.y0 < 0 or b.x1 > fig.bbox.width or b.y1 > fig.bbox.height:
        outside.append(t.get_text())
assert not outside, outside
assert not legend.get_window_extent(renderer).overlaps(
    ax.xaxis.label.get_window_extent(renderer)
)
assert selected[['ci95LowerInput', 'ci95LowerPM']].min().min() > -8
assert selected[['ci95UpperInput', 'ci95UpperPM']].max().max() < 8
stem = 'Figure_2C_top10_protein_coding_paired_effects_titled_3x4'
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
