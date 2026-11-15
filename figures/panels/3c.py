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

OUT = panel_output('3C')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '3C'

P = PROJECT_ROOT
source = P / 'results/transcriptFeatureAssociation.tsv'
digest = hashlib.sha256(source.read_bytes()).hexdigest()
table = pd.read_csv(source, sep = '\t')
order = [
    'logFivePrimeUtrLength',
    'fivePrimeUtrGc',
    'hasCandidateUorf',
    'hasOverlappingUorf',
    'kozakScore',
    'logCdsLength',
    'logThreePrimeUtrLength',
]
labels = [
    '5′ UTR length',
    '5′ UTR GC fraction',
    'Contained candidate\nuORF',
    'Overlapping candidate\nuORF',
    'Kozak score',
    'CDS length',
    '3′ UTR length',
]
d = (
    table[
        table.model.eq('joint')
        & table.analysis.eq('redistribution')
        & table.variant.eq('primary')
    ]
    .set_index('feature')
    .loc[order]
    .reset_index()
)
assert (
    len(d) == 7
    and d.adjustedPValue.lt(0.05).sum() == 5
    and d.transcripts.eq(11679).all()
    and d.genes.eq(3843).all()
)
assert (
    np.isfinite(
        d[
            [
                'standardizedCoefficient',
                'confidenceIntervalLow',
                'confidenceIntervalHigh',
                'adjustedPValue',
            ]
        ]
    )
    .all()
    .all()
)
ink = '#202938'
gray = '#63707D'
orange = '#D55E00'
blue = '#0072B2'
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
fig, ax = plt.subplots(figsize = (5, 6), dpi = 150)
fig.subplots_adjust(left = 0.35, right = 0.83, bottom = 0.22, top = 0.89)
fig.text(
    0.035, 0.963, 'C', ha = 'left', va = 'center', fontsize = 14, fontweight = 'bold'
)
fig.text(
    0.14,
    0.963,
    'mRNA architecture is associated with P/M differences',
    ha = 'left',
    va = 'center',
    fontsize = 11,
)
fig.text(0.975, 0.91, r'$P_{\mathrm{BH}}$', ha = 'right', va = 'center', fontsize = 9)
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = '#E7EBEF', linewidth = 0.55)
ax.axvline(0, color = gray, linewidth = 0.75, zorder = 1)
for i, r in d.iterrows():
    color = orange if r.standardizedCoefficient > 0 else blue
    ax.errorbar(
        r.standardizedCoefficient,
        i,
        xerr = [
            [r.standardizedCoefficient - r.confidenceIntervalLow],
            [r.confidenceIntervalHigh - r.standardizedCoefficient],
        ],
        fmt = 'o',
        color = color,
        ecolor = color,
        markerfacecolor = color if r.adjustedPValue < 0.05 else 'white',
        markeredgecolor = 'black',
        markeredgewidth = 0.6,
        markersize = 6,
        capsize = 2.5,
        elinewidth = 1.05,
        zorder = 3,
    )
    p = float(r.adjustedPValue)
    if p < 0.001:
        exponent = int(np.floor(np.log10(p)))
        mantissa = p / 10**exponent
        text = r'$' + f'{mantissa:.3g}' + r'\times10^{' + str(exponent) + r'}$'
    else:
        text = f'{p:.5f}'
    yy = 0.22 + 0.67 * (6.6 - i) / 7.2
    fig.text(0.975, yy, text, ha = 'right', va = 'center', fontsize = 8.5)
ax.set_ylim(6.6, -0.6)
ax.set_xlim(-0.12, 0.15)
ax.set_xticks([-0.10, 0, 0.10], ['−0.10', '0', '0.10'], fontsize = 8.5)
ax.set_yticks(range(7))
ax.tick_params(axis = 'y', length = 0, labelleft = False)
ax.spines['left'].set_visible(False)
label_transform = blended_transform_factory(fig.transFigure, ax.transData)
feature_labels = [
    ax.text(
        0.05,
        i,
        label,
        transform = label_transform,
        ha = 'left',
        va = 'center',
        fontsize = 9,
        clip_on = False,
    )
    for i, label in enumerate(labels)
]
ax.set_xlabel(r'Standardized $\beta$', labelpad = 8)
handles = [
    Line2D([], [], color = orange, linewidth = 1.4, label = 'Higher P/M in ASD'),
    Line2D([], [], color = blue, linewidth = 1.4, label = 'Lower P/M in ASD'),
    Line2D(
        [],
        [],
        marker = 'o',
        color = ink,
        linestyle = 'none',
        markeredgecolor = 'black',
        markersize = 5.5,
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
        markersize = 5.5,
        label = r'$P_{\mathrm{BH}}\geq0.05$',
    ),
]
legend = fig.legend(
    handles = handles,
    loc = 'lower center',
    bbox_to_anchor = (0.5, 0.018),
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
legend_box = legend.get_window_extent(renderer)
heading_x = legend_box.x0 / fig.bbox.width + 0.01
heading_y = legend_box.y1 / fig.bbox.height + 0.008
color_heading = fig.text(
    heading_x,
    heading_y,
    'Associated with:',
    ha = 'left',
    va = 'bottom',
    fontsize = 8.3,
)
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
assert not legend.get_window_extent(renderer).overlaps(
    ax.xaxis.label.get_window_extent(renderer)
)
assert all(
    t.get_window_extent(renderer).x1 < ax.bbox.x0 - 4 * fig.dpi / 72
    for t in feature_labels
)
assert (
    abs(
        (
            legend.get_window_extent(renderer).x0
            + legend.get_window_extent(renderer).x1
        )
        / 2
        - fig.bbox.width / 2
    )
    < 1
)
assert (
    d.confidenceIntervalLow.min() > ax.get_xlim()[0]
    and d.confidenceIntervalHigh.max() < ax.get_xlim()[1]
)
stem = 'Figure_3C_mRNA_architecture_standardized_beta_5x6'
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
