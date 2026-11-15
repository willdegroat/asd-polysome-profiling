import matplotlib
from _common import PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from reportlab.pdfgen import canvas

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('1legend')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '1legend'

plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'text.color': '#202938',
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    }
)
fig = plt.figure(figsize = (5, 0.4), dpi = 150)
fig.text(0.035, 0.75, 'Condition:', ha = 'left', va = 'center')
for x, c, label in [(0.24, '#238B45', 'Unaffected'), (0.50, '#C43C39', 'ASD')]:
    for offset, marker in [(0, 'o'), (0.028, '^'), (0.056, 's')]:
        fig.add_artist(
            Line2D(
                [x + offset],
                [0.75],
                transform = fig.transFigure,
                marker = marker,
                color = c,
                markeredgecolor = 'black',
                markeredgewidth = 0.5,
                markersize = 6,
                linestyle = 'none',
            )
        )
    fig.text(x + 0.086, 0.75, label, ha = 'left', va = 'center')
fig.text(0.035, 0.25, 'Library:', ha = 'left', va = 'center')
for x, m, label in [
    (0.24, 'o', 'Input'),
    (0.50, '^', 'Monosome'),
    (0.78, 's', 'Polysome'),
]:
    fig.add_artist(
        Line2D(
            [x],
            [0.25],
            transform = fig.transFigure,
            marker = m,
            color = '#63707D',
            markeredgecolor = 'black',
            markeredgewidth = 0.5,
            markersize = 6,
            linestyle = 'none',
        )
    )
    fig.text(x + 0.03, 0.25, label, ha = 'left', va = 'center')
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for t in fig.findobj(Text):
    if t.get_visible() and t.get_text().strip():
        b = t.get_window_extent(renderer)
        assert (
            b.x0 >= 0
            and b.y0 >= 0
            and b.x1 <= fig.bbox.width
            and b.y1 <= fig.bbox.height
        ), t.get_text()
stem = 'Figure_1_shared_PCA_legend_5x0p4'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 240)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((1200, 96))
    im.save(OUT / 'preview.png')
c = canvas.Canvas(str(OUT / f'{stem}.pdf'), pagesize = (360, 28.8))
c.drawImage(str(OUT / f'{stem}.png'), 0, 0, width = 360, height = 28.8)
c.showPage()
c.save()
