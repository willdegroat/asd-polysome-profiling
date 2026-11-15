import hashlib

import matplotlib
import numpy as np
import pandas as pd
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, to_rgb
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2G')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2G'

source = PROJECT_ROOT / 'results/fuAsdRankedEnrichment.tsv'
digest = hashlib.sha256(source.read_bytes()).hexdigest()
data = pd.read_csv(source, sep = '\t')
layers = ['input', 'monosome', 'polysome', 'redistribution']
tiers = ['Fu ASD FDR \u2264 0.001', 'Fu ASD FDR \u2264 0.05']
nes = data.pivot(
    index = 'tier', columns = 'context', values = 'normalizedEnrichmentScore'
).loc[tiers, layers]
q = data.pivot(
    index = 'tier', columns = 'context', values = 'globalAdjustedPValue'
).loc[tiers, layers]
assert nes.shape == (2, 4) and q.lt(0.05).to_numpy().sum() == 4
assert np.isfinite(nes.to_numpy()).all() and np.abs(nes.to_numpy()).max() <= 3
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
        'font.size': 9,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'text.color': ink,
        'axes.labelcolor': ink,
        'xtick.color': 'black',
        'ytick.color': 'black',
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (5, 1.5), dpi = 150)
fig.text(0.02, 0.98, 'G', ha = 'left', va = 'top', fontsize = 14, fontweight = 'bold')
fig.text(
    0.14,
    0.934,
    'ASD-associated genes favor monosome over polysome association in ASD',
    ha = 'left',
    va = 'center',
    fontsize = 8.75,
    fontweight = 'normal',
)
ax = fig.add_axes([0.15, 0.31, 0.84, 0.53])
ax.imshow(
    nes.to_numpy(), cmap = cmap, norm = norm, aspect = 'auto', interpolation = 'nearest'
)
ax.set_yticks([0, 1], ['ASD72', 'ASD185'], fontsize = 9)
ax.set_xticks(range(4), ['Input', 'Monosome', 'Polysome', 'P/M'], fontsize = 9)
ax.tick_params(length = 0, pad = 4)
for spine in ax.spines.values():
    spine.set_visible(False)


def luminance(rgb):
    a = np.asarray(rgb[:3])
    a = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
    return float(a @ np.array([0.2126, 0.7152, 0.0722]))


inklum = luminance(to_rgb(ink))
labels = []
for r in range(2):
    for c in range(4):
        value = float(nes.iloc[r, c])
        significant = bool(q.iloc[r, c] < 0.05)
        label = f'{value:+.2f}'.replace('-', '\u2212') + (
            '*' if significant else ''
        )
        lum = luminance(cmap(norm(value)))
        color = (
            'white'
            if 1.05 / (lum + 0.05) > (lum + 0.05) / (inklum + 0.05)
            else ink
        )
        ax.text(
            c,
            r,
            label,
            ha = 'center',
            va = 'center',
            fontsize = 10.4,
            fontweight = 'bold' if significant else 'normal',
            color = color,
        )
        labels.append(label)
fig.text(
    0.57,
    0.025,
    r'* $P_{\mathrm{BH}}<0.05$',
    ha = 'center',
    va = 'bottom',
    fontsize = 8,
    color = muted,
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
stem = 'Figure_2G_ASD_gene_enrichment_titled_5x1p5'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 900)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1500, 450))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 108))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 108)
c.showPage()
c.save()
