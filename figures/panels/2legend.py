import matplotlib
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('2legend')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '2legend'

palette = {'negative': '#526779', 'zero': '#F7F7F7', 'positive': '#B5367A'}
cmap = LinearSegmentedColormap.from_list(
    'NES_slate_white_magenta', list(palette.values())
)
norm = TwoSlopeNorm(vmin = -3, vcenter = 0, vmax = 3)
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 8,
        'text.color': 'black',
        'xtick.color': 'black',
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (5, 0.5), dpi = 150)
ax = fig.add_axes([0.07, 0.64, 0.86, 0.22])
bar = ColorbarBase(
    ax,
    cmap = cmap,
    norm = norm,
    orientation = 'horizontal',
    ticks = [-3, -2, -1, 0, 1, 2, 3],
)
bar.outline.set_edgecolor('black')
bar.outline.set_linewidth(0.5)
bar.ax.tick_params(labelsize = 8, length = 2, width = 0.5, pad = 1)
bar.ax.set_xticklabels(['\u22123', '\u22122', '\u22121', '0', '1', '2', '3'])
fig.text(
    0.5,
    0.045,
    'Normalized enrichment score (NES)',
    ha = 'center',
    va = 'bottom',
    fontsize = 8,
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
stem = 'Figure_2_GHI_shared_NES_bar_5x0p5'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 300)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1500, 150))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 36))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 36)
c.showPage()
c.save()
