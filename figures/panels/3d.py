import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('3D')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '3D'

d = pd.read_csv(HERE / 'Figure_3D_source_data.csv')
stats = pd.read_csv(
    PROJECT_ROOT / 'manuscript/followup/all_models_with_inference.csv'
)
r = stats[stats.model.eq('same_CDS_and_5UTR_sequence')].iloc[0]
x = d.withinMatchedGroupUTRlengthSD.to_numpy()
partial = d.adjustedPMpartialResidualLog2.to_numpy()
b = np.array([r.standardizedCoefficient])
sd = pd.Series({'log2FoldChange': r.outcomeSD})
blue = '#0072B2'
ink = '#202938'
gray = '#63707D'
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
fig.text(0.035, 0.955, 'D', fontsize = 14, fontweight = 'bold', va = 'center')
title = fig.text(
    0.14,
    0.955,
    'Longer 3′ UTRs are associated with lower P/M in ASD versus Unaffected',
    fontsize = 11,
    va = 'center',
)
fig.text(0.17, 0.885, 'Identical 5′ UTR and CDS', fontsize = 9, ha = 'left')
sk = fig.add_axes([0.17, 0.788, 0.78, 0.060])
sk.set_xlim(0, 1)
sk.set_ylim(0, 1)
sk.axis('off')
for y, end in [(0.68, 0.68), (0.20, 0.96)]:
    sk.add_patch(
        Rectangle(
            (0.03, y - 0.035), 0.16, 0.07, fc = 'white', ec = 'black', lw = 0.75
        )
    )
    sk.add_patch(
        Rectangle(
            (0.19, y - 0.095), 0.34, 0.19, fc = '#52606D', ec = 'black', lw = 0.75
        )
    )
    sk.add_patch(
        Rectangle(
            (0.53, y - 0.035), end - 0.53, 0.07, fc = 'white', ec = 'black', lw = 0.75
        )
    )
for xpos, label in [(0.11, '5′ UTR'), (0.36, 'CDS'), (0.74, '3′ UTR')]:
    sk.text(xpos, 1.02, label, ha = 'center', fontsize = 8, clip_on = False)
fig.text(0.17, 0.758, r'$P$ = 0.0118', fontsize = 8.5, color = blue)
ax = fig.add_axes([0.17, 0.095, 0.78, 0.645])
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = '#E7EBEF', lw = 0.55)
ax.axhline(0, color = '#C0C6CC', lw = 0.65)
ax.axvline(0, color = '#C0C6CC', lw = 0.65)
xx = np.linspace(min(x) - 0.12, max(x) + 0.12, 200)
slope = b[0] * sd.log2FoldChange
lo = r.confidenceIntervalLow * sd.log2FoldChange
hi = r.confidenceIntervalHigh * sd.log2FoldChange
ax.fill_between(
    xx,
    np.minimum(lo * xx, hi * xx),
    np.maximum(lo * xx, hi * xx),
    color = blue,
    alpha = 0.13,
    lw = 0,
    zorder = 1,
)
ax.scatter(
    x,
    partial,
    s = 17,
    facecolor = '#9AA6B0',
    edgecolor = 'black',
    linewidth = 0.4,
    alpha = 0.7,
    zorder = 2,
)
ax.plot(xx, slope * xx, color = blue, lw = 1.7, zorder = 3)
ax.set_xlim(xx.min(), xx.max())
pad = (max(partial) - min(partial)) * 0.07
ax.set_ylim(min(partial) - pad, max(partial) + pad)
ax.set_xlabel('Relative 3′ UTR length (SD)', labelpad = 7, fontsize = 9.5)
ax.set_ylabel(
    'Adjusted P/M change\n' + r'($\log_2$ ASD/Unaffected)',
    labelpad = 7,
    fontsize = 9.5,
)
ax.tick_params(labelsize = 8.5)
ax.set_xticks(
    np.arange(np.ceil(ax.get_xlim()[0]), np.floor(ax.get_xlim()[1]) + 1)
)
ax.set_yticks(
    np.arange(np.ceil(ax.get_ylim()[0]), np.floor(ax.get_ylim()[1]) + 1)
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
        box = t.get_window_extent(renderer)
        if (
            box.x0 < 0
            or box.y0 < 0
            or box.x1 > fig.bbox.width
            or box.y1 > fig.bbox.height
        ):
            outside.append(t.get_text())
assert not outside, outside
assert np.isfinite(partial).all()
assert ax.get_ylim()[0] < min(partial) and ax.get_ylim()[1] > max(partial)
stem = 'Figure_3D_matched_3UTR_association_clean_5x6'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 3600)
    im.save(OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600))
    im.thumbnail((850, 1020))
    im.save(OUT / 'preview.png')
