#!/usr/bin/env python3
"""
Render LAMMPS UFSW dump files with OVITO (no clicking needed).

Install once:   pip install ovito
Examples:
  python render_fsw.py                              # last frame, colour by type -> fsw.png
  python render_fsw.py --cu                         # plate B is Cu (Al/Cu run)
  python render_fsw.py --water hide                 # remove water
  python render_fsw.py --color temp                 # colour metal by temperature (v_Tat)
  python render_fsw.py --movie                      # every frame -> fsw.mp4
  python render_fsw.py --frame 0 --out start.png    # a specific frame
"""
import argparse
import numpy as np
from ovito.io import import_file
from ovito.modifiers import (SelectTypeModifier, DeleteSelectedModifier,
                             ComputePropertyModifier, ColorCodingModifier)
from ovito.vis import Viewport, TachyonRenderer

p = argparse.ArgumentParser()
p.add_argument("--dump", default="dump.fsw.*", help="dump file or wildcard pattern")
p.add_argument("--cu", action="store_true", help="type 4 is Cu (run made with -var cu 1)")
p.add_argument("--water", choices=["show", "transparent", "hide"], default="transparent")
p.add_argument("--color", choices=["type", "temp", "speed"], default="type")
p.add_argument("--tmax", type=float, default=900, help="upper end of temperature colour scale (K)")
p.add_argument("--frame", type=int, default=-1, help="frame to render (-1 = last)")
p.add_argument("--movie", action="store_true", help="render all frames to a movie")
p.add_argument("--out", default=None)
p.add_argument("--size", default="1600x900")
a = p.parse_args()

pipe = import_file(a.dump)
nframes = pipe.source.num_frames

# ---- names, radii, colours for the 4 atom types (done in code, not by hand) ----
TYPES = {1: ("Al plate A", 1.43, (0.66, 0.68, 0.71)),
         2: ("Tool (C)",   0.77, (0.15, 0.15, 0.15)),
         3: ("Water (mW)", 1.40, (0.37, 0.66, 0.91)),
         4: ("Cu plate B", 1.28, (0.85, 0.51, 0.23)) if a.cu else
            ("Al plate B", 1.43, (0.85, 0.51, 0.23))}

def style_types(frame, data):
    pt = data.particles_.particle_types_
    for tid, (name, r, c) in TYPES.items():
        t = pt.type_by_id_(tid)
        t.name, t.radius, t.color = name, r, c
pipe.modifiers.append(style_types)

# ---- water handling ----
if a.water == "hide":
    pipe.modifiers.append(SelectTypeModifier(property="Particle Type", types={3}))
    pipe.modifiers.append(DeleteSelectedModifier())
elif a.water == "transparent":
    pipe.modifiers.append(ComputePropertyModifier(
        output_property="Transparency", expressions=["0.8*(ParticleType==3)"]))

data0 = pipe.compute(0)

# ---- optional colour coding of the metal ----
if a.color != "type":
    names = {k.lower(): k for k in data0.particles.keys()}
    if a.color == "temp":
        if "v_tat" not in names:
            raise SystemExit("v_Tat not in dump - colour by temperature not available")
        prop, lo, hi = names["v_tat"], 250, a.tmax
    else:
        prop, lo, hi = "Velocity Magnitude", 0, 15     # A/ps
    pipe.modifiers.append(SelectTypeModifier(property="Particle Type", types={1, 4}))
    pipe.modifiers.append(ColorCodingModifier(property=prop, start_value=lo, end_value=hi,
                                              gradient=ColorCodingModifier.Jet(),
                                              only_selected=True))

pipe.add_to_scene()

# ---- camera ----
vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(1.0, 1.6, -0.9))
vp.zoom_all()

w, h = map(int, a.size.lower().split("x"))
ren = TachyonRenderer(ambient_occlusion=False, shadows=False)

if a.movie:
    out = a.out or "fsw.mp4"
    vp.render_anim(filename=out, size=(w, h), fps=10, renderer=ren)
else:
    frame = a.frame if a.frame >= 0 else nframes - 1
    out = a.out or "fsw.png"
    vp.render_image(filename=out, size=(w, h), frame=frame, renderer=ren)
print(f"wrote {out}  ({nframes} frames in dump)")
