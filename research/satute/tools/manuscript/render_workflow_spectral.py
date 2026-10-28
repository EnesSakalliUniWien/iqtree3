#!/usr/bin/env python3
"""Render a static model/spectral extension of Figure 1 and its offline HTML.

Run from research/satute: .venv/bin/python tools/manuscript/render_workflow_spectral.py
Reuses the approved workflow SVG, model helpers, and publication exporters.
"""
from __future__ import annotations

import html
import io
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import proj3d

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(PROJECT / 'src'))
from satute_analysis.models import build_q, reversible_eigendecomposition
from tools.figures.export import export_svg_pdf, export_svg_png
from tools.figures.latex_svg import compile_latex_svg, latex_svg_element, _prefix_ids
from tools.manuscript.render_workflow_latex import EQUATIONS

MODELS = (
    ('JC', 'JC', 'Equal exchangeabilities; equal frequencies'),
    ('K80', 'K2P{3}+FQ', 'Transitions α; transversions β'),
    ('GTR', 'GTR{1,4,0.7,1.2,3,0.8}+F{0.35,0.15,0.25,0.25}',
     'Unequal exchangeabilities and frequencies'),
)
STATES = ('A', 'C', 'G', 'T')
# Regular tetrahedron; probability vectors map affinely to its interior.
RADIUS = 1.45
BASE_RADIUS = RADIUS * np.sqrt(8 / 9)
VERTICES = np.array([[0, 0, RADIUS], [-BASE_RADIUS, 0, -RADIUS/3],
                     [BASE_RADIUS/2, BASE_RADIUS*np.sqrt(3)/2, -RADIUS/3],
                     [BASE_RADIUS/2, -BASE_RADIUS*np.sqrt(3)/2, -RADIUS/3]])
INK, MUTED, RUST = '#202020', '#56646c', '#8d2f14'
PATH_COLOUR = '#b55518'
MODE_COLOURS = ('#204e85', '#3d7778', '#775477')
ORANGE_INK = '#a44b0d'
GEOMETRY_INK = '#7c8992'
MODE_BAR_TIME = 1.0
WORKFLOW_MATH_SIZE = 8.9
ORIGINAL = PROJECT / 'artifacts/releases/generated-figures/figure1_satute_weighted_workflow.svg'
DEST = PROJECT / 'manuscript/figures'
STEM = 'figure1_satute_models_spectral'
HTML_PATH = DEST / 'figure1_satute_weighted_workflow.html'
PRINT_WIDTH_MM = 180
ANNOTATION_SIZE = 19
MODEL_WIDTH, MODEL_HEIGHT = 1296, 734.4
MODEL_LABEL_OFFSETS = (
    ((52, 52), (55, 26), (-48, 40)),
    ((-50, -32), (48, 40), (60, -28)),
    ((54, -34), (58, 30), (-55, 42)),
)
STATIONARY_LABEL_OFFSETS = ((-46, -36), (-45, 18), (-52, -34))

SYMBOLIC_MATRICES = (
    r'Q=\frac{1}{\mu}\begin{pmatrix}'
    r'-3\alpha&\alpha&\alpha&\alpha\\'
    r'\alpha&-3\alpha&\alpha&\alpha\\'
    r'\alpha&\alpha&-3\alpha&\alpha\\'
    r'\alpha&\alpha&\alpha&-3\alpha\end{pmatrix}',
    r'Q=\frac{1}{\mu}\begin{pmatrix}'
    r'-s&\beta&\alpha&\beta\\\beta&-s&\beta&\alpha\\'
    r'\alpha&\beta&-s&\beta\\\beta&\alpha&\beta&-s\end{pmatrix}',
    r'Q=\frac{1}{\mu}\begin{pmatrix}'
    r'-s_A&a\pi_C&b\pi_G&c\pi_T\\'
    r'a\pi_A&-s_C&d\pi_G&e\pi_T\\'
    r'b\pi_A&d\pi_C&-s_G&f\pi_T\\'
    r'c\pi_A&e\pi_C&f\pi_G&-s_T\end{pmatrix}',
)
SYMBOLIC_DEFINITIONS = (
    r'\mu=3\alpha',
    r's=\alpha+2\beta,\quad\mu=s',
    r'\textstyle s_i=\sum_{j\ne i}r_{ij}\pi_j,\quad\mu=\sum_i\pi_i s_i',
)


def math_element(expression, *, x, y, key, height=None, font_size=None, colour=INK):
    """Use the original workflow's LaTeX fonts and isolate repeated glyph IDs."""
    if font_size is not None:
        height = compile_latex_svg(expression, colour).height_pt * font_size / 10
    if height is None:
        raise ValueError('Supply a height or font_size for mathematical artwork')
    element = latex_svg_element(expression, x=x, y=y, height=height,
                                colour=colour, data_kind=key)
    return _prefix_ids(element, key)


def spectral_model(spec: str) -> dict:
    q, pi = build_q(spec)
    values, right, inverse = reversible_eigendecomposition(q, pi)
    order = np.argsort(values)[::-1]
    values = values[order]
    # Rows of inverse are left eigenvectors of Q; these become columns for Q^T.
    modes = inverse[order[1:], :].T.copy()
    for j in range(3):
        if modes[np.argmax(np.abs(modes[:, j])), j] < 0:
            modes[:, j] *= -1
    basis = np.column_stack((pi, modes))
    return {'q': q, 'pi': pi, 'values': values, 'modes': modes, 'basis': basis}


def probability(model: dict, p0: np.ndarray, t: float) -> np.ndarray:
    coefficients = np.linalg.solve(model['basis'], p0)
    decay = np.exp(model['values'] * t)
    decay[0] = 1.0
    return model['basis'] @ (coefficients * decay)


def mode_weights(model: dict, length: float) -> np.ndarray:
    values = model['values'][1:]
    return np.exp((values-values[0])*length)


MODE_BAR_WEIGHTS = mode_weights(spectral_model(MODELS[2][1]), MODE_BAR_TIME)


def plot_models() -> str:
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'DejaVu Sans'],
                         'font.size': 17, 'text.color': INK, 'mathtext.fontset': 'cm',
                         'svg.fonttype': 'path', 'svg.hashsalt': 'satute-spectral-figure'})
    fig = plt.figure(figsize=(18, 10.2), facecolor='white')
    overlays = []
    # One scale for all three symbolic matrices, including the wider GTR form.
    compiled_matrices = [compile_latex_svg(expr) for expr in SYMBOLIC_MATRICES]
    matrix_scale = min(2.0, 385 / max(m.width_pt for m in compiled_matrices))
    p0 = np.array([1., 0., 0., 0.])
    for col, (name, spec, subtitle) in enumerate(MODELS):
        model = spectral_model(spec)
        center = (col + .5) / 3
        overlays.append(f'<circle cx="{col*432+44}" cy="26" r="20.5" fill="#263238"/>')
        overlays.append(text(col*432+44,36.5,'abc'[col],31,'#ffffff','middle',True))
        fig.text((col*432+79)/MODEL_WIDTH, 1-35/MODEL_HEIGHT, name,
                 ha='left', fontsize=25, weight='bold')
        subtitle = subtitle.replace('; equal frequencies', ';\nequal frequencies').replace(
            'exchangeabilities and frequencies', 'exchangeabilities\nand frequencies')
        fig.text(center, 1-62/MODEL_HEIGHT, subtitle, ha='center', va='top',
                 fontsize=ANNOTATION_SIZE, color=MUTED, linespacing=1.15)
        overlays.append(math_element(SYMBOLIC_MATRICES[col],
                        x=center*MODEL_WIDTH, y=125, height=compiled_matrices[col].height_pt*matrix_scale,
                        key=f'matrix-{name}'))
        overlays.append(math_element(SYMBOLIC_DEFINITIONS[col], x=center*MODEL_WIDTH, y=238,
                        font_size=ANNOTATION_SIZE, key=f'normalization-{name}'))
        fig.text(center, 1-287/MODEL_HEIGHT, 'Exact JC spectrum' if name=='JC' else 'Illustrative parameter values',
                 ha='center', fontsize=ANNOTATION_SIZE, color=MUTED)
        if name=='JC':
            pi_expression = r'\pi=(\tfrac14,\tfrac14,\tfrac14,\tfrac14)^{\mathsf T}'
            eigen_expression = r'\lambda_1=\lambda_2=\lambda_3=-\frac{4\alpha}{\mu}=-\frac43'
        else:
            pi = ', '.join(f'{x:.2f}' for x in model['pi'])
            values = ', '.join(f'{x:.3f}' for x in model['values'][1:])
            pi_expression = r'\pi = ('+pi+r')^{\mathsf{T}}'
            eigen_expression = r'(\lambda_1,\lambda_2,\lambda_3)=('+values+r')'
        overlays.append(math_element(pi_expression, x=center*MODEL_WIDTH, y=304,
                        font_size=ANNOTATION_SIZE, key=f'pi-{name}'))
        overlays.append(math_element(eigen_expression, x=center*MODEL_WIDTH, y=340,
                        font_size=ANNOTATION_SIZE, key=f'eigenvalues-{name}'))
        repeated = ['One decay rate; multiplicity 3', 'Two decay rates; multiplicities 1 + 2',
                    'Three distinct decay rates'][col]
        fig.text(center, 1-397/MODEL_HEIGHT, repeated,
                 ha='center', va='top', fontsize=ANNOTATION_SIZE, color=MUTED, linespacing=1.1)
        ax = fig.add_axes([col/3+.005, .045, .323, .43], projection='3d', computed_zorder=False)
        ax.set_proj_type('ortho'); ax.view_init(elev=17, azim=-69)
        ax.set_box_aspect((1,1,1), zoom=1.24); ax.set_axis_off()
        for i in range(4):
            for j in range(i+1,4):
                ax.plot(*VERTICES[[i,j]].T, color=GEOMETRY_INK, lw=1.3, zorder=1)
        for i, vertex in enumerate(VERTICES):
            ax.scatter(*vertex, s=23, color=INK, depthshade=False, zorder=5)
            offset=np.array([0,0,.14]) if i==0 else (vertex * .12)
            ax.text(*(vertex+offset), STATES[i], fontsize=19, weight='bold', ha='center', zorder=10)
        stationary = model['pi'] @ VERTICES
        label_points = []
        for j, colour in enumerate(MODE_COLOURS):
            v = model['modes'][:,j]
            # Both endpoints stay inside the simplex; lengths are display-scaled.
            scale = .88 * np.min(model['pi'][np.abs(v)>1e-12]/np.abs(v[np.abs(v)>1e-12]))
            delta = (scale*v) @ VERTICES
            ax.plot(*np.vstack((stationary-delta,stationary+delta)).T,
                    color=colour, ls='--', lw=1.6, zorder=4)
            ax.quiver(*stationary,*delta, color=colour, linewidth=1.7, arrow_length_ratio=.2, zorder=5)
            label_points.append(stationary+delta)
        points = np.array([probability(model,p0,t) @ VERTICES for t in np.linspace(0,12,240)])
        ax.plot(*points.T, color=PATH_COLOUR, linewidth=3, zorder=6)
        for time,marker in ((.5,'o'),(2.,'s')):
            point=probability(model,p0,time) @ VERTICES
            ax.scatter(*point,s=40,marker=marker,color=PATH_COLOUR,edgecolors='white',linewidth=.6,
                       depthshade=False,zorder=7)
        ax.scatter(*stationary,s=75,marker='D',color=INK,edgecolors='white',depthshade=False,zorder=8)
        ax.set(xlim=(-1.55,1.1),ylim=(-1.45,1.45),zlim=(-.85,1.75))
        for j, (point, offset, colour) in enumerate(zip(label_points, MODEL_LABEL_OFFSETS[col], MODE_COLOURS)):
            px,py,_=proj3d.proj_transform(*point,ax.get_proj())
            ax.annotate(f'$v_{j+1}$', (px,py), xytext=offset, textcoords='offset points',
                        color=colour, fontsize=ANNOTATION_SIZE+1, ha='center', va='center',
                        arrowprops={'arrowstyle':'-', 'color':colour, 'lw':.7}, zorder=12)
        px,py,_=proj3d.proj_transform(*stationary,ax.get_proj())
        ax.annotate(r'$\pi$', (px,py), xytext=STATIONARY_LABEL_OFFSETS[col],textcoords='offset points',
                    color=INK,fontsize=22,ha='center',va='center',
                    arrowprops={'arrowstyle':'-', 'color':INK, 'lw':.8},zorder=12)
        note=['Straight trajectory\nto the stationary distribution', 'Curved trajectory\nwith two decay rates',
              'Curved trajectory\nwith three decay rates'][col]
        fig.text(center,.012,note,ha='center',va='bottom',fontsize=ANNOTATION_SIZE,
                 color=MUTED,linespacing=1.1)
    buffer=io.StringIO(); fig.savefig(buffer,format='svg',metadata={'Date':None}); plt.close(fig)
    raw=buffer.getvalue()
    return raw[raw.index('<svg'):].replace('</svg>', '\n'.join(overlays)+'</svg>', 1)


def nested(svg: str, x: float, y: float, width: float, height: float, view_box: str | None=None) -> str:
    opening=svg[:svg.index('>')+1]
    opening=re.sub(r'\s(?:width|height|x|y)=[\"\'][^\"\']*[\"\']','',opening)
    if view_box:
        opening=re.sub(r'viewBox=[\"\'][^\"\']*[\"\']',f'viewBox="{view_box}"',opening)
    opening=opening.replace('<svg ',f'<svg x="{x}" y="{y}" width="{width}" height="{height}" overflow="hidden" ',1)
    return opening+svg[svg.index('>')+1:]


def text(x: float,y: float,content: str,size=24,colour=INK,anchor='start',bold=False) -> str:
    return (f'<text x="{x}" y="{y}" font-family="Helvetica, Arial, sans-serif" font-size="{size}" '
            f'fill="{colour}" text-anchor="{anchor}" font-weight="{"bold" if bold else "normal"}">'
            f'{html.escape(content)}</text>')


def refine_workflow_text(svg: str) -> str:
    """Adjust output typography and labels without changing the source equations."""
    labels = {
        'ML-tree inference': 'Maximum-likelihood tree inference',
        'branch-length weights': 'relative persistence',
        'site-pattern coefficient': 'weighted coherence',
    }

    def refine_label(match):
        element = match.group(0)
        content = html.unescape(re.sub(r'<[^>]*>', '', element))
        coordinates = {name: float(value) for name, value in
                       re.findall(r'\b(x|y)=[\"\']([\d.-]+)[\"\']', element)}
        x, y = coordinates.get('x', -1), coordinates.get('y', -1)
        if content in {'mode 1', 'mode 2', 'mode 3'}:
            return ''
        if content in {'relative persistence', 'over AB'} and x > 450:
            return ''
        # Enlarge complete likelihood tables together, preserving row spacing.
        for center, left, right in ((328.19, 314, 343), (379.82, 365, 395)):
            if left < x < right and 121 < y < 148:
                if 'Courier' in element:
                    element = re.sub(r'\bx=[\"\'][\d.-]+[\"\']', f'x="{x-1.2:.2f}"', element, count=1)
                return (f'<g transform="translate({center} 121.46) scale(1.65) '
                        f'translate({-center} -121.46)">{element}</g>')
        if content in labels:
            element = element[:element.index('>')+1] + html.escape(labels[content]) + '</text>'
        if 'Helvetica' in element:
            element = re.sub(r'font-size:\s*([\d.]+)px',
                             lambda m: f'font-size: {max(8.8, float(m[1])):.2f}px', element)
            # Let glyphs use their actual metrics after changing size or wording.
            element = re.sub(r'\s(?:textLength|lengthAdjust)=[\"\'][^\"\']*[\"\']', '', element)
        return element

    svg = re.sub(r'<text\b[^>]*>.*?</text>', refine_label, svg, flags=re.DOTALL)

    def remove_schematic_bar(match):
        element = match.group(0)
        points = re.search(r'points=[\"\']([^\"\']+)', element)
        if points:
            coordinates = [tuple(map(float,pair.split(','))) for pair in points[1].split()]
            if all(x>419 and 348<y<385 for x,y in coordinates):
                return ''
        return element
    svg = re.sub(r'<polygon\b[^>]*/>', remove_schematic_bar, svg)

    def enlarge_table(match):
        element = match.group(0)
        center = 328.19 if "x='314.13'" in element else 379.82
        return (f'<g transform="translate({center} 121.46) scale(1.65) '
                f'translate({-center} -121.46)">{element}</g>')

    svg = re.sub(r"<rect x='(?:314\.13|365\.76)'[^>]*/>", enlarge_table, svg)
    # Make room beside the enlarged left likelihood table.
    def move_subtree_label(match):
        return re.sub(r'\bx="([\d.-]+)"', lambda m: f'x="{float(m[1])-18:.6f}"',
                      match.group(0), count=1)
    svg = re.sub(r'<svg\b[^>]*data-kind="latex-left-subtree"[^>]*>', move_subtree_label, svg)
    return svg.replace('#777777', MUTED).replace('#DC840C', ORANGE_INK)


def enlarge_workflow_math(svg: str) -> str:
    sizes = {spec.kind: spec.font_size for spec in EQUATIONS}

    def resize(match):
        opening = match.group(0)
        kind = re.search(r'data-kind="([^"]+)"', opening)[1]
        if kind not in sizes:
            return opening
        values = {key: float(value) for key,value in
                  re.findall(r'\b(x|y|width|height)="([\d.-]+)"', opening)}
        scale = max(1, WORKFLOW_MATH_SIZE/sizes[kind])
        cx,cy = values['x']+values['width']/2, values['y']+values['height']/2
        cy += {'latex-weighted-coherence': 3, 'latex-alignment-average': 3}.get(kind, 0)
        if kind=='latex-critical-value':
            cy=452
        values.update(x=cx-values['width']*scale/2, y=cy-values['height']*scale/2,
                      width=values['width']*scale, height=values['height']*scale)
        return re.sub(r'\b(x|y|width|height)="[\d.-]+"',
                      lambda m:f'{m[1]}="{values[m[1]]:.6f}"', opening)

    svg=re.sub(r'<svg\b[^>]*data-kind="[^"]+"[^>]*>', resize, svg)
    # Side-set summation limits preserve the equations while leaving breathing room.
    for spec in EQUATIONS:
        if spec.kind not in {'latex-weighted-coherence','latex-alignment-average'}:
            continue
        expression=r'\textstyle '+spec.expression
        height=compile_latex_svg(expression,spec.colour).height_pt*WORKFLOW_MATH_SIZE/10
        replacement=math_element(expression,x=spec.x*511,
                                 y=(1-spec.y)*504+3-height/2,height=height,key=spec.kind,
                                 colour=spec.colour)
        pattern=rf'<svg\b[^>]*data-kind="{re.escape(spec.kind)}"[^>]*>.*?</svg>'
        svg=re.sub(pattern,lambda _:replacement,svg,flags=re.DOTALL)
    return svg


def weight_bars() -> str:
    parts=[f'<g role="img" aria-label="Illustrated GTR relative persistence at t=1" '
           f'data-model="{html.escape(MODELS[2][1], quote=True)}" data-branch-length="{MODE_BAR_TIME}">',
           text(445,324,'Illustrated GTR',8.8,MUTED,'middle'),
           text(445,337,'t = 1',8.8,MUTED,'middle')]
    for i,(weight,colour) in enumerate(zip(MODE_BAR_WEIGHTS,MODE_COLOURS)):
        baseline=354+15*i
        parts += [text(417,baseline,f'mode {i+1}',8.8,MUTED,'end'),
                  f'<rect x="421" y="{baseline-5.5}" width="{42*weight:.8f}" height="5" '
                  f'style="fill:{colour};stroke:none" data-weight="{weight:.12f}" data-mode="{i+1}"/>',
                  text(485,baseline,f'{weight:.2f}',8.8,INK,'end')]
    return '\n'.join(parts+['</g>'])


def layout_workflow(svg: str) -> str:
    """Recompose the four existing panels on a shared print grid.

    Artwork and equations retain their source coordinates and a uniform scale;
    the panel frames and workflow connectors have independent reserved space.
    """
    root = ET.fromstring(svg)
    artwork = root.find(".//*[@class='svglite']")
    drawing = list(artwork)[-1]
    children = list(drawing)
    starts = [i for i, node in enumerate(children[:-1])
              if node.tag.endswith('circle') and children[i+1].text in ('d', 'e', 'f', 'g')]
    if len(starts) != 4:
        raise ValueError('Expected the four established workflow panels')
    panels = {}
    for index, start in enumerate(starts):
        stop = starts[index+1] if index+1 < len(starts) else len(children)
        panels[children[start+1].text] = children[start+3:stop]
    equations = {node.get('data-kind'): node for node in root[-1]}
    kinds = {
        'd': ('latex-site-pattern',),
        'e': ('latex-left-subpattern', 'latex-right-subpattern', 'latex-left-subtree',
              'latex-right-subtree', 'latex-left-likelihood', 'latex-right-likelihood',
              'latex-coherence-coefficient'),
        'f': ('latex-persistence-weight', 'latex-weighted-coherence', 'latex-alignment-average',
              'latex-weighted-statistic', 'latex-null-distribution'),
        'g': ('latex-critical-value', 'latex-alpha', 'latex-decision-rule'),
    }
    # 76.4 mm frames, a 10 mm column gutter, and an 8 mm connector lane.
    # Uniform 3x artwork scale raises the base mathematics to 7.57 pt in print.
    frames = (
        ('d', 'Alignment and site pattern', 86, 1490, 570, 14.56, 61.74),
        ('e', 'Branch likelihood at AB', 950, 1490, 570, 237, 65.74),
        ('g', 'Decision under the null', 86, 2140, 560, 14.56, 290),
        ('f', 'Weighted coherence and z-score', 950, 2140, 560, 251, 290),
    )
    parts = ['<g data-kind="workflow-print-grid">']
    # Keep the established SVG drawing rules, but not the old clipping frame.
    parts.extend(ET.tostring(node, encoding='unicode') for node in artwork
                 if node.tag.endswith('style') or node.tag.endswith('defs'))
    scale, border_width, bottom_inset = 3, 2.25, 30  # 3 mm at 180 mm print width.
    bottom_equations = {'e': 'latex-coherence-coefficient',
                        'f': 'latex-null-distribution', 'g': 'latex-decision-rule'}
    for key, title, x, y, height, source_x, source_y in frames:
        # Anchor the equation viewport to one shared INNER-border inset.
        bottom_shift = 0
        if key in bottom_equations:
            equation = equations[bottom_equations[key]]
            target_bottom = source_y + (height-border_width/2-bottom_inset)/scale
            bottom_shift = target_bottom-float(equation.get('y'))-float(equation.get('height'))
        parts += [f'<rect id="workflow-frame-{key}" x="{x}" y="{y}" width="764" height="{height}" rx="8" '
                  f'fill="#FBFCFD" stroke="#D7DEE8" stroke-width="{border_width}"/>',
                  f'<circle cx="{x+43}" cy="{y+48}" r="29" fill="#263238"/>',
                  text(x+43, y+59, key, 43, '#ffffff', 'middle', True),
                  text(x+91, y+57, title, 36, INK, bold=True),
                  f'<g class="svglite" transform="translate({x} {y}) scale(3) translate({-source_x} {-source_y})">']
        def place(element: str, dy: float = 0) -> str:
            if key == 'd':
                # Enlarge the sparse alignment uniformly, centred below its heading.
                return ('<g transform="translate(141.893333 174) scale(1.12) '
                        f'translate(-137 -150.5)">{element}</g>')
            return f'<g transform="translate(0 {dy:.8f})">{element}</g>'

        def branch_y(value: float) -> float:
            # Reserve a table-to-branch gap while retaining all glyph sizes.
            return (value-189)*0.52+177.4

        for index, node in enumerate(panels[key]):
            dy = bottom_shift
            if key == 'e':
                dy = -15
                if index < 25:
                    # Change drawing coordinates, never scale text or circle radii.
                    node = ET.fromstring(ET.tostring(node, encoding='unicode'))
                    for attr in ('y', 'y1', 'y2', 'cy'):
                        if node.get(attr) is not None:
                            node.set(attr, str(branch_y(float(node.get(attr)))))
                    if node.get('points'):
                        node.set('points', ' '.join(
                            f'{px},{branch_y(float(py)):.8f}' for px, py in
                            (pair.split(',') for pair in node.get('points').split())))
                    if node.text in ('A', 'B'):
                        node.set('y', '167.4')
                    elif node.text == 'AB':
                        node.set('y', '192.4')
                    dy = 0
                if node.text == 'coherence coefficients':
                    dy = -10 + bottom_shift + 6
            parts.append(place(ET.tostring(node, encoding='unicode'), dy))
        for kind in kinds[key]:
            dy = bottom_shift
            if key == 'e' and kind != 'latex-coherence-coefficient':
                if kind in ('latex-left-likelihood', 'latex-right-likelihood'):
                    dy = -15
                elif kind in ('latex-left-subpattern', 'latex-right-subpattern'):
                    dy = branch_y(155.99)-3-float(equations[kind].get('height'))-float(equations[kind].get('y'))
                else:
                    dy = branch_y(144.14)-15-float(equations[kind].get('y'))
            parts.append(place(ET.tostring(equations[kind], encoding='unicode'), dy))
        if key == 'f':
            parts.append(place(weight_bars(), bottom_shift))
        parts.append('</g>')
    # Connectors live in explicit lanes and terminate at their destination frames.
    parts += ['<defs><marker id="workflow-arrow" markerWidth="9" markerHeight="9" '
              'refX="8" refY="4.5" orient="auto" markerUnits="userSpaceOnUse">'
              '<path d="M0,0 L8,4.5 L0,9" fill="white" stroke="#53636f" stroke-width="1.5"/>'
              '</marker></defs>',
              '<g fill="none" stroke="#53636f" stroke-width="2.4" marker-end="url(#workflow-arrow)">'
              '<path d="M800,1490 V1473 H1000 V1490"/>'
              '<path d="M1675,2060 V2138"/>'
              '<path d="M950,2630 H852"/>'
              '<path d="M125,2140 V2062"/></g>',
              text(900,1453,'Maximum-likelihood tree inference',27,MUTED,'middle'),
              text(1645,2109,'SatuTe',29,RUST,'end',True),
              text(160,2109,'data re-evaluation',29,RUST,bold=True),
              '</g>']
    return '\n'.join(parts)


def outline_text(svg: str) -> str:
    """Use the installed exporter so the distributable SVG needs no fonts."""
    executable=shutil.which('inkscape')
    if executable is None:
        raise RuntimeError('Portable SVG export requires the installed Inkscape CLI')
    # The R drawing has no native paths in this group. Once text becomes paths,
    # its generic path rule would otherwise apply black strokes to every glyph.
    svg=svg.replace('.svglite path,', '')
    with tempfile.TemporaryDirectory(prefix='satute-outline-') as folder:
        source=Path(folder)/'source.svg'; destination=Path(folder)/'outlined.svg'
        source.write_text(svg)
        subprocess.run([executable,str(source),'--export-text-to-path','--export-plain-svg',
                        f'--export-filename={destination}'],check=True,capture_output=True,text=True)
        result=destination.read_text()
    root=ET.fromstring(result)
    if any(e.tag.rsplit('}',1)[-1] in {'text','tspan'} for e in root.iter()):
        raise RuntimeError('Portable SVG still contains font-dependent text')
    return result[result.index('<svg'):].strip()


CAPTION_TITLE = 'Figure 1. Substitution models, spectral decay and SatuTe.'
CAPTION_PARAGRAPHS = (
    ('a–c. Substitution models.',
     'The symbolic rate matrices use state order A, C, G, T. Each row sums to zero; μ normalizes '
     'the stationary mean substitution rate to one. In JC, α is the common rate before normalization. '
     'In K80, α and β are transition and transversion rates before normalization, with s = α + 2β. '
     'In GTR, (a, b, c, d, e, f) denote the symmetric exchangeabilities '
     '(r_AC, r_AG, r_AT, r_CG, r_CT, r_GT), and s_i = Σ_{j≠i} r_ij π_j.'),
    ('Illustrative parameter values.',
     'JC frequencies and eigenvalues are exact. K80 uses α/β = 3. GTR uses exchangeabilities '
     '(1, 4, 0.7, 1.2, 3, 0.8) in the order above and stationary frequencies '
     '(0.35, 0.15, 0.25, 0.25). The numerical K80 and GTR spectra and trajectories illustrate '
     'these choices; they are not fitted manuscript estimates.'),
    ('3D probability space.',
     'Each tetrahedron is a 3D affine representation of four-state probability space. '
     'Probability vectors are columns: dp/dt = Qᵀp and Qᵀπ = 0. '
     'The nonstationary eigenvectors satisfy Qᵀvᵢ = λᵢvᵢ. The orange trajectory starts at '
     'p(0) = (1, 0, 0, 0)ᵀ. Circles mark t = 0.5, squares mark t = 2, and the black diamond '
     'marks the stationary distribution π as t → ∞. Time is measured in expected substitutions '
     'per site at stationarity. Dashed eigenvector directions are translated to π and scaled '
     'for display; they are not probability trajectories. Eigenvector signs and bases within '
     'repeated eigenspaces are arbitrary, while the probability trajectory is independent of '
     'that choice. JC has a multiplicity-three nonstationary eigenvalue; K80 here has multiplicities '
     'one and two.'),
    ('d–g. SatuTe.',
     'The manuscript uses h_i for left eigenvectors of Q. With the same eigenbasis and '
     'normalization, h_i = v_i: the probability directions above are the coherence directions below. '
     'They are normalized by h_iᵀ diag(π)⁻¹ h_j = δ_ij. '
     '(d) Alignment columns define site patterns. (e) Branch AB divides each site pattern '
     'into two subpatterns; their partial likelihoods give coherence coefficients. '
     '(f) Relative persistence weights compare each component with the slowest nonstationary '
     'component. SatuTe uses the weighted coherence to construct a z-score. '
     f'The mode bars use the illustrated GTR model at t = {MODE_BAR_TIME:g}; '
     f'the relative weights are ({", ".join(f"{w:.3f}" for w in MODE_BAR_WEIGHTS)}). '
     '(g) Rejection of subtree independence identifies phylogenetic signal.'),
)
CAPTION = CAPTION_TITLE + ' ' + ' '.join(f'{label} {body}' for label, body in CAPTION_PARAGRAPHS)


def render() -> None:
    model_svg=plot_models()
    original=ORIGINAL.read_text(); original=original[original.index('<svg'):].strip()
    original=re.sub(r'>([abcd])</text>',lambda m:'>'+dict(zip('abcd','defg'))[m[1]]+'</text>',original)
    original=enlarge_workflow_math(original)
    original=refine_workflow_text(original)
    workflow=layout_workflow(original)
    parts=['<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
           f'width="{PRINT_WIDTH_MM}mm" height="272mm" viewBox="0 0 1800 2720" role="img" aria-labelledby="title" aria-describedby="description">',
           '<title id="title">Substitution models, spectral decay and SatuTe</title>',
           '<desc id="description">'+html.escape(CAPTION)+'</desc>',
           '<rect width="1800" height="2720" fill="white"/>',
           text(900,57,'Substitution models, spectral decay and SatuTe',43,RUST,'middle',True),
           text(900,103,'Rate matrices (A, C, G, T), spectral components, probability trajectories and SatuTe',25,MUTED,'middle'),
           nested(model_svg,20,120,1760,1760*MODEL_HEIGHT/MODEL_WIDTH),
           math_element(r'Q^{\mathsf T}=V\,\operatorname{diag}(0,\lambda_1,\lambda_2,\lambda_3)V^{-1},'
                        r'\qquad V=[\,\pi\ v_1\ v_2\ v_3\,]',
                        x=900,y=1127,font_size=27,key='spectral-decomposition'),
           math_element(r'\textstyle p(t)=\pi+\sum_{i=1}^{3}c_i e^{\lambda_i t}v_i,'
                        r'\qquad\lambda_0=0,\quad\lambda_1,\lambda_2,\lambda_3<0',
                        x=900,y=1170,font_size=27,key='probability-expansion'),
           '<line x1="260" y1="1230" x2="320" y2="1230" stroke="'+PATH_COLOUR+'" stroke-width="4"/>',
           text(334,1238,'p(t), initial state A',26),
           '<circle cx="635" cy="1230" r="6" fill="'+PATH_COLOUR+'"/>',text(651,1238,'t = 0.5',26),
           '<rect x="850" y="1224" width="12" height="12" fill="'+PATH_COLOUR+'"/>',text(878,1238,'t = 2',26),
           '<path d="M 1090,1230 l 8,-8 l 8,8 l -8,8 Z" fill="'+INK+'"/>',text(1122,1238,'stationary distribution π',26),
           text(900,1280,'Dashed arrows show eigenvector directions at π; lengths are scaled for display.',26,MUTED,'middle'),
           '<line x1="86" y1="1310" x2="1714" y2="1310" stroke="#bdc7cd" stroke-width="1.5"/>',
           text(900,1343,'SatuTe uses the same left eigenvectors and relative persistence',26,RUST,'middle'),
           math_element(r'v_i=h_i,\qquad Q^{\mathsf T}h_i=\lambda_i h_i',
                        x=900,y=1355,font_size=27,key='probability-coherence-directions'),
           math_element(r'w_i(t)=\exp\{(\lambda_i-\lambda_1)t\},\qquad w_1(t)=1',
                        x=900,y=1392,font_size=27,key='relative-persistence',colour=RUST),
           workflow,
           '</svg>']
    svg=outline_text('\n'.join(parts))
    DEST.mkdir(parents=True,exist_ok=True)
    svg_path=DEST/f'{STEM}.svg'; svg_path.write_text(svg+'\n')
    export_svg_pdf(svg_path,DEST/f'{STEM}.pdf')
    export_svg_png(svg_path,DEST/f'{STEM}.png',width=1800)
    page='''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Figure 1 — Substitution models, spectral decay and SatuTe</title>
<style>
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { margin: 0; background: #f2f3f3; color: #202020; font-family: Helvetica, Arial, sans-serif; }
main { max-width: 1400px; margin: auto; padding: 24px; }
figure { margin: 0; }
figure > svg { display: block; width: 100%; height: auto; background: white; }
figcaption { max-width: 70ch; margin: 28px auto; font-size: 16px; line-height: 1.55; }
figcaption p { margin: 0 0 1em; }
::selection { background: #dce8ee; color: #202020; }
@media (max-width: 600px) { main { padding: 8px; } figcaption { margin: 16px 8px; } }
@media print { @page { size: A4; margin: 10mm; } body { background: white; }
 main { max-width: none; padding: 0; } figure > svg { width: 180mm; max-width: 100%; height: auto; }
 figcaption { break-before: page; font-size: 10pt; line-height: 1.45; } }
</style></head><body><main><figure>
'''+svg+'\n<figcaption><p><strong>'+html.escape(CAPTION_TITLE)+'</strong></p>' + ''.join(
        '<p><strong>'+html.escape(label)+'</strong> '+html.escape(body)+'</p>'
        for label, body in CAPTION_PARAGRAPHS)+'</figcaption>\n</figure></main></body></html>\n'
    HTML_PATH.write_text(page)
    for path in [HTML_PATH,svg_path,DEST/f'{STEM}.pdf',DEST/f'{STEM}.png']:
        print(path)


if __name__=='__main__':
    render()
