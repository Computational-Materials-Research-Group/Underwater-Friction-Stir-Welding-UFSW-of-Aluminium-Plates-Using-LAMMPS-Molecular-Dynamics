#!/usr/bin/env python3
"""
Build ready-made OVITO session files (.ovito) for the UFSW dumps, so the OVITO
desktop app opens with names, radii, colours, transparency and render settings
already set. Optional colourings are included as modifiers that are switched OFF -
tick their checkbox in the pipeline to use them.

Run it in the folder with the dump files, using the SAME OVITO version as your
desktop app (check Help > About in OVITO, then:  pip install ovito==3.x.y):
    python3 make_ovito_session.py           # Al/Al run
    python3 make_ovito_session.py --cu      # Al/Cu run
Writes:
    fsw_atoms.ovito     atoms as spheres
    fsw_surface.ovito   solid surfaces (plates, tool, water) instead of atoms
Open with: OVITO > File > Load Session State.
"""
import argparse, glob, os
from ovito.io import import_file
from ovito.modifiers import (SelectTypeModifier, DeleteSelectedModifier, SliceModifier, ComputePropertyModifier, ColorCodingModifier,
                             ConstructSurfaceModifier, ClearSelectionModifier)
from ovito.vis import Viewport, TachyonRenderer
import ovito

p = argparse.ArgumentParser()
p.add_argument("--dump", default="dump.fsw.*")
p.add_argument("--cu", action="store_true", help="type 4 is Cu")
p.add_argument("--tmax", type=float, default=900)
a = p.parse_args()
if not glob.glob(a.dump):
    raise SystemExit(f"no files match {a.dump} - run this in the folder with the dump files")
dump = os.path.abspath(a.dump)

TYPES = {1: ("Al plate A", 1.43, (0.66, 0.68, 0.71)),
         2: ("Tool (C)",   0.77, (0.20, 0.20, 0.22)),
         3: ("Water (mW)", 1.40, (0.37, 0.66, 0.91)),
         4: ("Cu plate B", 1.28, (0.85, 0.51, 0.23)) if a.cu else
            ("Al plate B", 1.43, (0.85, 0.51, 0.23))}

def titled(mod, title, enabled=True):
    try: mod.title = title
    except Exception: pass
    mod.enabled = enabled
    return mod

def new_pipeline():
    for pl in list(ovito.scene.pipelines):
        pl.remove_from_scene()
    pipe = import_file(dump)
    # names / radii / colours stored in the file source -> no Python needed in the GUI
    for tid, (name, r, c) in TYPES.items():
        t = pipe.source.data.particles_.particle_types_.type_by_id_(tid)
        t.name, t.radius, t.color = name, r, c
    return pipe

def tool_center(pipe):
    d = pipe.compute(0)
    pos = d.particles.positions[d.particles.particle_types[...] == 2]
    return pos[:, 0].mean(), pos[:, 1].mean()

def render_setup(fname, cam=(1.0, 1.6, -0.9)):
    try:   # make the active (big) viewport a perspective view looking at the tool
        vc = ovito.scene.viewport_config
        vp = vc.active_vp or vc.viewports[-1]
        vp.type = Viewport.Type.Perspective
        vp.camera_dir = cam
        vp.zoom_all()
        vc.maximized_vp = vp
    except Exception as e:
        print("  (viewport layout not set:", e, ")")
    rs = ovito.scene.render_settings
    rs.renderer = TachyonRenderer(ambient_occlusion=True, shadows=False)
    rs.size = (1920, 1080)
    rs.output_filename = os.path.abspath(fname)
    rs.save_to_file = True
    try: ovito.scene.anim.current_frame = ovito.scene.anim.last_frame   # open on the final frame
    except Exception: pass

# ------------------------------------------------------------ 1) atoms session
pipe = new_pipeline()
xt, yt = tool_center(pipe)
m = pipe.modifiers
# water cutaway: remove the FRONT half of the water only (metal untouched), keep the back half
m.append(titled(SelectTypeModifier(property="Particle Type", types={3}), "Water cutaway (1/2): select water"))
m.append(titled(SliceModifier(normal=(0, -1, 0), distance=-yt, only_selected=True),
                "Water cutaway (2/2): remove front half of water"))
m.append(titled(ClearSelectionModifier(), "clear selection"))
m.append(titled(ComputePropertyModifier(output_property="Transparency",
                                        expressions=["0.35*(ParticleType==3)"]),
                "Water transparency (0 = solid)"))
m.append(titled(SelectTypeModifier(property="Particle Type", types={3}), "OPTIONAL: hide all water (1/2)", False))
m.append(titled(DeleteSelectedModifier(), "OPTIONAL: hide all water (2/2)", False))
names = {k.lower(): k for k in pipe.compute(0).particles.keys()}
m.append(titled(SelectTypeModifier(property="Particle Type", types={1, 4}), "OPTIONAL: select metal (for colouring)", False))
if "v_tat" in names:
    m.append(titled(ColorCodingModifier(property=names["v_tat"], start_value=250, end_value=a.tmax,
                                        gradient=ColorCodingModifier.Jet(), only_selected=True),
                    "OPTIONAL: colour metal by temperature", False))
m.append(titled(ColorCodingModifier(property="Velocity Magnitude", start_value=0, end_value=15,
                                    gradient=ColorCodingModifier.Jet(), only_selected=True),
                "OPTIONAL: colour metal by speed", False))
m.append(titled(ClearSelectionModifier(), "clear selection (keeps atoms from turning red)"))
pipe.add_to_scene()
render_setup("fsw_atoms_render.png", cam=(0.35, 1.0, -0.6))   # look from the cut-away side
ovito.scene.save("fsw_atoms.ovito")
print("wrote fsw_atoms.ovito")

# ------------------------------------------------------------ 2) surface session
pipe = new_pipeline()
m = pipe.modifiers
for tid in (1, 4, 2, 3):
    name, r, col = TYPES[tid]
    m.append(titled(SelectTypeModifier(property="Particle Type", types={tid}), f"select {name}"))
    cs = ConstructSurfaceModifier(method=ConstructSurfaceModifier.Method.AlphaShape,
                                  radius=4.0 if tid != 2 else 2.5, smoothing_level=8, only_selected=True,
                                  identify_regions=False)
    cs.vis.surface_color = col
    cs.vis.show_cap = True            # closes the surface where it meets the box edge -> solid look
    cs.vis.surface_transparency = 0.6 if tid == 3 else 0.0
    cs.vis.cap_color = col                         # filled cross-sections in the material's colour
    cs.vis.cap_transparency = 0.6 if tid == 3 else 0.0
    m.append(titled(cs, f"surface: {name}"))
m.append(titled(ClearSelectionModifier(), "clear selection"))
pipe.compute(0).particles.vis.enabled = False
pipe.source.data.particles.vis.enabled = False       # hide atoms, show only surfaces
pipe.add_to_scene()
render_setup("fsw_surface_render.png")
ovito.scene.save("fsw_surface.ovito")
print("wrote fsw_surface.ovito")
