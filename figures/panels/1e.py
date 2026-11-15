import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('1E')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '1E'

scores = pd.read_csv(HERE / 'PCA_scores.csv', dtype = {'sampleId': str})
coordinates = scores[['PC1', 'PC2']].to_numpy()
variance = np.array([43.2649202351, 15.6511147594])
ink = '#202938'
grid = '#E7EBEF'
colors = {'Unaffected': '#238B45', 'ASD': '#C43C39'}
markers = {'Input': 'o', 'Monosome': '^', 'Polysome': 's'}
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.edgecolor': 'black',
        'axes.linewidth': 0.75,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': ink,
        'axes.labelcolor': ink,
        'xtick.color': 'black',
        'ytick.color': 'black',
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig, ax = plt.subplots(figsize = (5, 2.9), dpi = 150)
fig.subplots_adjust(left = 0.14, right = 0.97, bottom = 0.20, top = 0.86)
fig.text(
    0.035, 0.95, 'E', ha = 'left', va = 'center', fontsize = 14, fontweight = 'bold'
)
fig.text(
    0.14,
    0.95,
    'Principal component analysis of all RNA fractions',
    ha = 'left',
    va = 'center',
    fontsize = 11,
    fontweight = 'normal',
)
ax.set_axisbelow(True)
ax.grid(axis = 'y', color = grid, linewidth = 0.6)
for row in scores.itertuples():
    ax.scatter(
        row.PC1,
        row.PC2,
        color = colors[row.condition],
        marker = markers[row.fraction],
        s = 55,
        edgecolors = 'black',
        linewidths = 0.6,
        zorder = 3,
    )
ax.margins(x = 0.20, y = 0.24)
ax.set_xlabel(f'PC1 ({variance[0]:.1f}%)', labelpad = 5)
ax.set_ylabel(f'PC2 ({variance[1]:.1f}%)', labelpad = 6)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
placed = []
labels = []
point_pixels = ax.transData.transform(coordinates)
radius = np.sqrt(55) / 2 * fig.dpi / 72 + 0.6 * fig.dpi / 72
for row in scores.itertuples():
    for dx, dy in [
        (6, 7),
        (6, -10),
        (-6, 7),
        (-6, -10),
        (10, 0),
        (-10, 0),
        (0, 12),
        (0, -14),
    ]:
        text = ax.annotate(
            row.displayLabel,
            (row.PC1, row.PC2),
            xytext = (dx, dy),
            textcoords = 'offset points',
            ha = 'left' if dx > 0 else 'right' if dx < 0 else 'center',
            va = 'center',
            fontsize = 9,
            zorder = 4,
        )
        fig.canvas.draw()
        box = text.get_window_extent(renderer).expanded(1.06, 1.12)
        point_collision = (
            (point_pixels[:, 0] + radius > box.x0)
            & (point_pixels[:, 0] - radius < box.x1)
            & (point_pixels[:, 1] + radius > box.y0)
            & (point_pixels[:, 1] - radius < box.y1)
        ).any()
        inside = (
            box.x0 > ax.bbox.x0
            and box.x1 < ax.bbox.x1
            and box.y0 > ax.bbox.y0
            and box.y1 < ax.bbox.y1
        )
        if (
            inside
            and not point_collision
            and not any(box.overlaps(b) for b in placed)
        ):
            placed.append(box)
            labels.append(text)
            break
        text.remove()
    else:
        raise AssertionError(
            'No collision-free placement for ' + row.displayLabel
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
stem = 'Figure_1E_all_fractions_PCA_5x2p9'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 1740)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1200, 720))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 208.8))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 208.8)
c.showPage()
c.save()
