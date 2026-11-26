#!/usr/bin/env python3
"""
Build the phylograftR hex logo (SVG master + PNG).

Usage (from the package root):
    python data-raw/make_logo.py --font path/to/Montserrat[wght].ttf

Outputs (default): man/figures/logo.svg and man/figures/logo.png

Requirements:
    pip install fonttools            # always needed (text is converted to outlines)
    pip install cairosvg             # optional, only to also write the PNG
    Montserrat variable font (OFL licence): https://fonts.google.com/specimen/Montserrat
    The font is only needed to rebuild the logo; the SVG contains no live text.

Design (matches the deepSTRAPP hex style):
    navy #002060 outline, tree and text; blue #0B64B5 for the "R"; red #FF0000 for
    grafted taxa and dashed time rings; white -> #2A7AC6 vertical gradient background.
    Circular time-calibrated fan phylogeny, grafted clades and isolated grafted tips
    in red, root stem reaching the hex tip, package name curved over the fan.
"""
import argparse
import math
import os
import sys

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont
from fontTools.varLib import instancer

# ---- palette ---------------------------------------------------------------
NAVY, BLUE, RED, BOTTOM = "#002060", "#0B64B5", "#FF0000", "#2A7AC6"

# ---- layout (user units, 200 x 200 working canvas) --------------------------
HEX = "100,3 183.6,51.5 183.6,148.5 100,197 16.4,148.5 16.4,51.5"
VIEWBOX = "12 -1.5 176 203"          # tight crop around the hexagon
CX, CY, R = 100, 117, 60             # fan centre and radius
LINE_W, GAP_DEG = 3.1, 24            # branch width; gap (deg) left open at the bottom
TEXT_RADIUS, TEXT_WIDTH, TEXT_TRACK = R + 14, 142, 40
RING_FRACTIONS = (0.36, 0.72)        # dashed "time slice" rings, as fractions of R


class Node:
    def __init__(self, h, *children, grafted=False):
        self.h, self.ch, self.g = h, list(children), grafted


def tip(grafted=False):
    return Node(0, grafted=grafted)


def isolated_graft(h):
    """A new tip grafted as sister of an existing tip."""
    return Node(h, tip(), tip(True))


def build_tree():
    """Ultrametric tree (h = node height, 1 = root). Grafted tips/clades are red.
    Node heights are chosen so that arcs avoid the dashed rings (h = 0.64 and 0.28)."""
    G = isolated_graft
    return Node(
        1.0,
        Node(.58,
             Node(.43, Node(.14, tip(), G(.07)), tip()),
             Node(.46, G(.1), Node(.2, tip(), Node(.09, tip(True), tip(True), grafted=True)))),
        Node(.72,
             Node(.40, tip(), Node(.16, tip(), tip())),
             Node(.5,
                  Node(.19, G(.08), tip()),
                  Node(.39, tip(), Node(.12, G(.06), Node(.06, tip(True), tip(True), grafted=True))))))


def fan_svg(cx=CX, cy=CY, radius=R, lw=LINE_W, gap=GAP_DEG):
    tree = build_tree()
    tips = []

    def collect(n):
        if not n.ch:
            tips.append(n)
        for c in n.ch:
            collect(c)

    collect(tree)
    span = 360 - gap
    for i, n in enumerate(tips):
        n.a = math.radians(90 + gap / 2 + span * i / (len(tips) - 1))

    def layout(n):
        for c in n.ch:
            layout(c)
        if n.ch:
            n.a = sum(c.a for c in n.ch) / len(n.ch)
        n.r = radius * (1 - n.h)

    layout(tree)
    P = lambda r, a: (cx + r * math.cos(a), cy + r * math.sin(a))
    out = []
    for fr in RING_FRACTIONS:
        out.append(f'<circle cx="{cx}" cy="{cy}" r="{radius * fr:.1f}" fill="none" stroke="{RED}" '
                   f'stroke-width="1.3" stroke-dasharray="3.2 3.2" stroke-linecap="butt"/>')

    def draw(n):
        for c in n.ch:
            col = RED if c.g else NAVY
            x1, y1 = P(n.r, n.a)
            xa, ya = P(n.r, c.a)
            if n.r > 0.5 and abs(c.a - n.a) > 1e-6:
                sweep = 1 if c.a > n.a else 0
                out.append(f'<path d="M{x1:.2f} {y1:.2f} A{n.r:.2f} {n.r:.2f} 0 0 {sweep} {xa:.2f} {ya:.2f}" '
                           f'fill="none" stroke="{col}" stroke-width="{lw}"/>')
            xc, yc = P(c.r, c.a)
            out.append(f'<line x1="{xa:.2f}" y1="{ya:.2f}" x2="{xc:.2f}" y2="{yc:.2f}" '
                       f'stroke="{col}" stroke-width="{lw}"/>')
            draw(c)
        if not n.ch:
            x, y = P(n.r, n.a)
            out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.5" fill="{RED if n.g else NAVY}"/>')

    draw(tree)
    # root stem down to the hex tip: a "grounded" tree
    out.append(f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="196" stroke="{NAVY}" '
               f'stroke-width="{lw * 1.25}" stroke-linecap="butt"/>')
    return "\n".join(out)


def arc_text_svg(glyphs, cmap, text, cx, cy, rb, width, colors, track, center_deg=-90):
    """Text converted to outlines and laid along a circular arc (letters upright, outward)."""
    adv = [glyphs[cmap[ord(c)]].width for c in text]
    sc = width / sum(adv)
    total = (sum(adv) + track * (len(text) - 1)) * sc
    out, pos = [], 0.0
    for i, ch in enumerate(text):
        mid = (pos + adv[i] / 2 * sc) - total / 2
        th = math.radians(center_deg) + mid / rb
        ph = th + math.pi / 2
        px, py = cx + rb * math.cos(th), cy + rb * math.sin(th)
        a, b = math.cos(ph) * sc, math.sin(ph) * sc
        e = px - math.cos(ph) * sc * adv[i] / 2
        f = py - math.sin(ph) * sc * adv[i] / 2
        pen = SVGPathPen(glyphs)
        glyphs[cmap[ord(ch)]].draw(TransformPen(pen, (a, b, b, -a, e, f)))
        out.append(f'<path d="{pen.getCommands()}" fill="{colors.get(i, NAVY)}"/>')
        pos += (adv[i] + track) * sc
    return "\n".join(out)


def build_svg(font_path):
    font = instancer.instantiateVariableFont(TTFont(font_path), {"wght": 800})
    glyphs, cmap = font.getGlyphSet(), font.getBestCmap()
    name = "phylograftR"
    text = arc_text_svg(glyphs, cmap, name, CX, CY, TEXT_RADIUS, TEXT_WIDTH,
                        {len(name) - 1: BLUE}, TEXT_TRACK)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="{VIEWBOX}">
<title>phylograftR hex logo</title>
<defs><linearGradient id="bg" x1="0" y1="3" x2="0" y2="197" gradientUnits="userSpaceOnUse">
<stop offset="0" stop-color="#FFFFFF"/><stop offset="0.30" stop-color="#FFFFFF"/><stop offset="0.9" stop-color="{BOTTOM}"/></linearGradient></defs>
<polygon points="{HEX}" fill="url(#bg)" stroke="{NAVY}" stroke-width="5.5" stroke-linejoin="miter"/>
<g stroke-linecap="round" stroke-linejoin="round">
{text}
{fan_svg()}
</g>
</svg>
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--font", required=True, help="Path to Montserrat[wght].ttf (variable font)")
    ap.add_argument("--outdir", default=os.path.join("man", "figures"), help="Output directory")
    ap.add_argument("--png-width", type=int, default=700, help="PNG width in pixels")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    svg_path = os.path.join(args.outdir, "logo.svg")
    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write(build_svg(args.font))
    print("wrote", svg_path)

    try:
        import cairosvg
        png_path = os.path.join(args.outdir, "logo.png")
        cairosvg.svg2png(url=svg_path, write_to=png_path, output_width=args.png_width)
        print("wrote", png_path)
    except ImportError:
        print("cairosvg not installed: PNG skipped (pip install cairosvg, or export logo.svg yourself)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
