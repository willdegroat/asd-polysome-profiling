import json
import re

import matplotlib
import numpy as np
from _common import ASSETS, PROJECT_ROOT, SOURCE_DATA, panel_output
from matplotlib.patches import (
    Circle,
    Ellipse,
    FancyBboxPatch,
    PathPatch,
    Rectangle,
)
from matplotlib.path import Path as MPath
from matplotlib.text import Text
from PIL import Image

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = panel_output('1A')
ROOT = OUT
REV = OUT
P = PROJECT_ROOT
PROJECT = PROJECT_ROOT
HERE = SOURCE_DATA / '1A'

ink = '#202938'
gray = '#414C58'
line = '#7B8791'
light = '#F3F5F7'
green = '#238B45'
red = '#C43C39'
navy = '#073763'
neutral = '#405C70'
plt.rcParams.update(
    {
        'font.family': 'Arial',
        'font.size': 9,
        'text.color': ink,
        'svg.fonttype': 'none',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
    }
)
fig = plt.figure(figsize = (5, 12), dpi = 150, facecolor = 'white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 5)
ax.set_ylim(12, 0)
ax.axis('off')


def text(x, y, s, size = 9, ha = 'center', color = ink, weight = 'normal', **kw):
    return ax.text(
        x,
        y,
        s,
        fontsize = size,
        ha = ha,
        va = 'center',
        color = color,
        fontweight = weight,
        linespacing = 1.2,
        **kw,
    )


BOX_BOUNDS = []
ARROW_SEGMENTS = []


def box(x, y, w, h, fc = light, ec = line, rounding = 0.055, lw = 0.7):
    BOX_BOUNDS.append((x - 0.025, y - 0.025, x + w + 0.025, y + h + 0.025))
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle = f'round,pad=0.025,rounding_size={rounding}',
        facecolor = fc,
        edgecolor = ec,
        linewidth = lw,
    )
    ax.add_patch(patch)
    return patch


def arrow(x, y, xx, yy, color = gray):
    ARROW_SEGMENTS.append(((x, y), (xx, yy)))
    return ax.annotate(
        '',
        xy = (xx, yy),
        xytext = (x, y),
        arrowprops = {
            'arrowstyle': '-|>',
            'color': color,
            'lw': 0.8,
            'shrinkA': 1,
            'shrinkB': 1,
            'mutation_scale': 8,
        },
    )


def head(y, n, label):
    box(0.22, y, 4.56, 0.38, navy, navy)
    text(0.40, y + 0.19, n, 10, ha = 'left', color = 'white', weight = 'bold')
    text(0.68, y + 0.19, label, 10, ha = 'left', color = 'white')


def marker(x, y, m, c, s = 31):
    ax.scatter(
        [x],
        [y],
        marker = m,
        c = c,
        s = s,
        edgecolors = 'black',
        linewidths = 0.55,
        zorder = 5,
    )


source = json.loads(
    (ASSETS / 'source-icon-geometry.json').read_text(encoding = 'utf-8-sig')
)
native = {r['shape']['shapeType']: r for r in source['icons']}
hc = native['ELLIPSE']
bc = native['ROUND_RECTANGLE']
native_width = hc['size']['width']['magnitude'] * hc['transform']['scaleX']
native_gap = bc['transform']['translateY'] - (
    hc['transform']['translateY'] + native_width
)
gap_ratio = native_gap / native_width
assert 0 < gap_ratio < 1


def person(x, y, color):
    side = 0.22
    radius = side / 2
    ax.add_patch(Circle((x, y), radius, fc = color, ec = 'none'))
    body_y = y + radius + side * gap_ratio
    ax.add_patch(
        FancyBboxPatch(
            (x - radius, body_y),
            side,
            side,
            boxstyle = 'round,pad=0,rounding_size=.035',
            fc = color,
            ec = 'none',
        )
    )


def ipsc_colony(x, y):
    ax.add_patch(
        Ellipse((x, y + 0.12), 0.98, 0.38, fc = 'white', ec = line, lw = 0.8)
    )
    for dx, dy in [
        (-0.22, 0.03),
        (-0.07, -0.08),
        (0.10, -0.08),
        (0.25, 0.02),
        (0.10, 0.15),
        (-0.08, 0.15),
        (0.02, 0.04),
    ]:
        ax.add_patch(
            Circle((x + dx, y + dy), 0.105, fc = '#D0DCE6', ec = neutral, lw = 0.65)
        )
        ax.add_patch(
            Circle(
                (x + dx + 0.012, y + dy - 0.012), 0.035, fc = '#748C9E', ec = 'none'
            )
        )


def curve(points, lw = 1):
    path = MPath(points, [MPath.MOVETO] + [MPath.CURVE4] * (len(points) - 1))
    ax.add_patch(
        PathPatch(path, fill = False, color = neutral, lw = lw, capstyle = 'round')
    )


def neuron(x, y):
    # Soma, nucleus, branched dendrites, and a long axon with terminals.
    ax.add_patch(
        Ellipse((x, y), 0.22, 0.18, angle = -18, fc = '#B8CCD9', ec = neutral, lw = 0.9)
    )
    ax.add_patch(Circle((x - 0.015, y - 0.008), 0.037, fc = neutral, ec = 'none'))
    for pts in [
        [(-0.07, -0.05), (-0.20, -0.15), (-0.17, -0.28), (-0.29, -0.35)],
        [(0.02, -0.08), (0.06, -0.20), (0.17, -0.23), (0.17, -0.39)],
        [(0.09, -0.02), (0.23, -0.08), (0.28, -0.18), (0.42, -0.22)],
        [(-0.10, 0.02), (-0.25, 0.04), (-0.28, 0.15), (-0.43, 0.16)],
    ]:
        curve([(x + dx, y + dy) for dx, dy in pts], 1.05)
    for pts in [
        [(-0.18, -0.21), (-0.28, -0.24), (-0.36, -0.20), (-0.43, -0.25)],
        [(-0.24, -0.30), (-0.21, -0.35), (-0.22, -0.42), (-0.18, -0.45)],
        [(0.12, -0.28), (0.23, -0.30), (0.28, -0.37), (0.35, -0.37)],
        [(0.30, -0.16), (0.37, -0.11), (0.43, -0.11), (0.48, -0.13)],
        [(-0.29, 0.12), (-0.35, 0.20), (-0.32, 0.25), (-0.38, 0.29)],
    ]:
        curve([(x + dx, y + dy) for dx, dy in pts], 0.75)
    curve(
        [
            (x + 0.03, y + 0.07),
            (x + 0.08, y + 0.34),
            (x + 0.42, y + 0.18),
            (x + 0.53, y + 0.34),
        ],
        1.05,
    )
    for pts in [
        [(0.50, 0.30), (0.59, 0.31), (0.63, 0.29), (0.69, 0.32)],
        [(0.53, 0.34), (0.53, 0.39), (0.59, 0.43), (0.64, 0.43)],
        [(0.40, 0.23), (0.46, 0.20), (0.52, 0.22), (0.57, 0.19)],
    ]:
        curve([(x + dx, y + dy) for dx, dy in pts], 0.75)
    for dx, dy in [(0.69, 0.32), (0.64, 0.43), (0.57, 0.19)]:
        ax.add_patch(Circle((x + dx, y + dy), 0.018, fc = neutral, ec = 'none'))


text(0.175, 0.28, 'A', 14, ha = 'left', weight = 'bold')
title = text(2.50, 0.28, 'Study design and analysis overview', 11)
text(2.5, 0.70, 'One sibling pair', 10)
person(1.15, 1.02, green)
person(3.85, 1.02, red)
ax.plot([1.40, 3.60], [1.18, 1.18], color = line, lw = 0.8)
text(1.15, 1.60, 'Unaffected', 9.5, color = green)
text(3.85, 1.60, 'ASD', 9.5, color = red)
for x in [1.15, 3.85]:
    arrow(x, 1.75, x, 1.93)
ax.plot([1.15, 3.85], [1.96, 1.96], color = line, lw = 0.8)
arrow(2.5, 1.96, 2.5, 2.13)
ipsc_colony(1.43, 2.48)
neuron(3.44, 2.49)
arrow(2.02, 2.47, 2.90, 2.47)
text(2.46, 2.24, 'NGN2 induction', 8)
text(1.43, 3.06, 'iPSCs', 9.5)
text(3.57, 3.06, 'Induced neurons', 9.5)
text(2.5, 3.31, '3 independently processed cultures per individual', 8.5)
arrow(2.5, 3.46, 2.5, 3.61)
text(2.5, 3.77, 'Matched RNA pools from each culture', 9.5)
for x, m, label, definition in [
    (1.03, 'o', 'Input', 'Unfractionated'),
    (2.50, '^', 'Monosome', 'Pooled 40Sâ€“80S'),
    (3.97, 's', 'Polysome', 'Pooled polysomes'),
]:
    for yy, c in [(4.02, green), (4.23, red)]:
        for dx in [-0.27, 0, 0.27]:
            marker(x + dx, yy, m, c)
    text(x, 4.48, label, 9)
    text(x, 4.71, definition, 8, color = gray)
ax.plot([1.03, 1.03, 3.97, 3.97], [4.87, 4.97, 4.97, 4.87], color = line, lw = 0.7)
ax.plot([2.50, 2.50], [4.87, 4.97], color = line, lw = 0.7)
arrow(2.5, 4.97, 2.5, 5.10)
text(2.5, 5.25, '18 paired-end RNA-seq libraries', 10)
text(2.5, 5.47, '2 individuals Ã— 3 cultures Ã— 3 fractions', 8.5, color = gray)
ax.plot([0.23, 4.77], [5.66, 5.66], color = line, lw = 0.6)
text(2.5, 5.88, 'RNA-seq processing and quantification', 10)
text(2.5, 6.10, 'nf-core/rnaseq', 8, color = gray)
for x, label, software in [
    (1.0, 'Read QC', 'FastQC / MultiQC'),
    (2.5, 'Alignment', 'STAR Â· GRCh38'),
    (4.0, 'Quantification', 'Salmon Â· GENCODE v50'),
]:
    text(x, 6.33, label, 9)
    text(x, 6.55, software, 7.5, color = gray)
arrow(1.42, 6.33, 1.88, 6.33)
arrow(3.02, 6.33, 3.39, 6.33)
arrow(2.5, 6.70, 2.5, 6.85)
for r in range(6):
    for c in range(18):
        shade = ['#D3DDE5', '#8AA0B0', '#405F78', '#B4C3CF'][
            (r * 3 + c * 2) % 4
        ]
        ax.add_patch(
            Rectangle(
                (0.62 + c * 0.08, 6.94 + r * 0.073),
                0.071,
                0.063,
                fc = shade,
                ec = 'white',
                lw = 0.2,
            )
        )
text(1.35, 7.53, 'Gene / transcript counts', 8.5)
text(2.50, 7.15, '+', 13, color = gray)
box(3.00, 6.945, 1.30, 0.42, '#E5EDF3', line)
for r, c in enumerate([green, red]):
    yy = 7.04 + r * 0.21
    ax.plot([3.08, 3.35], [yy, yy], color = c, lw = 2.5)
    for xx, m in [(3.60, 'o'), (3.88, '^'), (4.15, 's')]:
        marker(xx, yy, m, c, 16)
text(3.65, 7.53, 'Sample metadata', 8.5)
text(2.5, 7.82, 'Expression filtering Â· library QC Â· PCA', 8.5)
ax.plot([1.00, 3.35], [7.94, 7.94], color = line, lw = 0.75)
ax.plot([2.5, 2.5], [7.91, 7.94], color = line, lw = 0.75)
for x in [1.00, 3.35]:
    arrow(x, 7.94, x, 8.09)
box(0.30, 8.16, 1.40, 0.35, '#4B5661', '#4B5661', rounding = 0.025)
text(1.00, 8.335, 'Abundance', 10, color = 'white')
text(1.00, 8.69, 'PyDESeq2', 9, color = gray)
text(1.00, 8.93, 'Condition effects\nwithin each pool', 8.5)
arrow(1.0, 9.20, 1.0, 9.52)
text(1.00, 9.68, 'Input RNA abundance', 8.5)
arrow(1.0, 9.84, 1.0, 10.12)
text(1.00, 10.30, 'GO overrepresentation', 8)
arrow(1.0, 10.45, 1.0, 10.61)
pts = np.array(
    [[0.65, 10.20], [0.94, 10.06], [1.31, 10.16], [0.78, 10.48], [1.18, 10.46]]
) + np.array([0, 0.70])
for a, b in [(0, 1), (1, 2), (0, 3), (1, 3), (1, 4), (2, 4), (3, 4)]:
    ax.plot(pts[[a, b], 0], pts[[a, b], 1], color = line, lw = 0.7)
for j, (x, y) in enumerate(pts):
    ax.add_patch(
        Circle(
            (x, y), 0.045 if j == 1 else 0.035, fc = neutral, ec = 'white', lw = 0.4
        )
    )
text(1.00, 11.42, 'Neural pathways', 9)
box(1.94, 8.16, 2.82, 0.35, '#234E64', '#234E64', rounding = 0.025)
text(3.35, 8.335, 'Redistribution', 10, color = 'white')
text(3.35, 8.69, 'Monopoly', 10, color = ink)
text(3.35, 8.95, 'Paired negative-binomial model', 9)
text(3.35, 9.17, 'Culture + fraction + condition Ã— fraction', 7.8)
text(3.35, 9.40, 'Fraction-wise normalization', 8, color = gray)
box(1.94, 9.53, 2.82, 0.35, '#263544', '#263544', rounding = 0.02, lw = 0.65)
text(
    3.35,
    9.705,
    "model.contrast('poly', 'mono')",
    8.5,
    color = 'white',
    fontfamily = 'Consolas',
)
text(3.35, 10.06, 'ASDâ€“Unaffected change in P/M', 8.5)
eq = text(
    3.35,
    10.42,
    r'$\log_2\left[\frac{(P/M)_{\mathrm{ASD}}}{(P/M)_{\mathrm{Unaffected}}}\right]$',
    16,
    math_fontfamily = 'dejavusans',
)
text(3.35, 10.80, r'Effects Â· uncertainty Â· adjusted $P$ values', 8)
text(
    3.35,
    11.00,
    'Sensitivity: culture omission and normalization',
    7.5,
    color = gray,
)
ax.plot(
    [2.47, 2.47, 4.23, 4.23], [11.21, 11.12, 11.12, 11.21], color = line, lw = 0.7
)
ax.plot([3.35, 3.35], [11.07, 11.12], color = line, lw = 0.7)
for x in [2.47, 4.23]:
    arrow(x, 11.21, x, 11.29)
for x, label, detail, ref in [
    (2.47, 'Gene-level', 'SynGO / ASD genes', 'Figure 2Câ€“I'),
    (4.23, 'Transcript-level', 'UTRs / CDS / uORFs', 'Figure 3'),
]:
    text(x, 11.42, label, 9)
    text(x, 11.65, detail, 7.8, color = gray)
for start, end in ARROW_SEGMENTS:
    for xmin, ymin, xmax, ymax in BOX_BOUNDS:
        for x, y in np.linspace(start, end, 101):
            assert not (xmin <= x <= xmax and ymin <= y <= ymax), (
                start,
                end,
                (xmin, ymin, xmax, ymax),
            )
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
outside = []
labels = [t.get_text() for t in fig.findobj(Text) if t.get_visible()]
assert "model.contrast('poly', 'mono')" in labels
assert '2 individuals Ã— 3 cultures Ã— 3 fractions' in labels
assert 'nf-core/rnaseq' in labels
assert not any(
    'Individuals' in label or 'Cultures' in label or 'Fractions' in label
    for label in labels
)
assert all(
    t.get_ha() == 'center'
    for t in fig.findobj(Text)
    if t.get_visible() and t.get_text().strip() and t.get_text() != 'A'
)
assert not any(
    re.search(r'Figure [23]|PCA \(B', t.get_text()) for t in fig.findobj(Text)
)
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
for label in fig.findobj(Text):
    if label.get_text() in ['PyDESeq2', 'Monopoly']:
        assert (
            label.get_window_extent(renderer).y1
            < (12 - 8.535) / 12 * fig.bbox.height - 1
        )
stem = 'Figure_1A_study_design_analysis_polished_5x12'
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi = 600)
plt.close(fig)
with Image.open(OUT / f'{stem}.png') as im:
    assert im.size == (3000, 7200)
    im.convert('RGB').save(
        OUT / f'{stem}.tiff', compression = 'tiff_lzw', dpi = (600, 600)
    )
    im.thumbnail((750, 1800))
    im.save(OUT / 'preview.png')
