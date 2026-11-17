# Underwater Friction Stir Welding (UFSW) of Aluminium Plates Using LAMMPS Molecular Dynamics

![LAMMPS](https://img.shields.io/badge/LAMMPS-Molecular%20Dynamics-blue?style=for-the-badge)
![Process](https://img.shields.io/badge/Process-Underwater%20FSW-purple?style=for-the-badge)
![Metal](https://img.shields.io/badge/pair__style-eam%2Falloy-orange?style=for-the-badge)
![Water](https://img.shields.io/badge/Water-mW%20(Stillinger--Weber)-teal?style=for-the-badge)
![OVITO](https://img.shields.io/badge/OVITO-Visualization-green?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24296019-blue?style=for-the-badge)](https://doi.org/10.5281/zenodo.24296019)

An **atomistic simulation of underwater friction stir welding (UFSW)** of two butted aluminium plates (or a dissimilar Al/Cu pair), using LAMMPS. A rigid diamond tool with a pin and a shoulder plunges into the seam between the plates, dwells, and then traverses along the seam while rotating. The plates are covered by a layer of coarse-grained mW water, which carries heat away from the weld zone. The same input runs without water (`-var wet 0`) as the dry reference.

The simulation captures **frictional and plastic heating under the tool**, **stirring and mixing of material across the weld seam**, **heat transfer into the surrounding water**, **tool forces and torque in each phase**, and **the crystal structure (defects and recrystallisation) of the stirred zone**. Temperature, forces, torque and the number of mixed atoms are written every 1 ps, and the full trajectory is written for visualisation in OVITO.

<img width="1280" height="720" alt="fsw_surface_render" src="https://github.com/user-attachments/assets/6413fc85-999a-4d15-8031-4cf05edb29ec" />


---

## Physics

The simulation models the weld as a molecular dynamics system with four atom types:

- **Plates:** two fcc plates (Al/Al, or Al/Cu with `-var cu 1`) butted along a seam in the x direction, with a 2 Å clearance at the seam
- **Metal interactions:** all metal–metal forces (Al–Al, Al–Cu, Cu–Cu) come from an embedded-atom (EAM) potential
- **Tool:** a rigid body of diamond-lattice carbon atoms whose motion (rotation + plunge + traverse) is prescribed by `fix move`; tool–tool forces are switched off
- **Tool–metal contact:** Morse potentials, so friction and plastic work come from atomic interactions, not from a friction law
- **Water:** the mW model, one bead per molecule, with a three-body Stillinger–Weber term that favours tetrahedral (hydrogen-bond-like) angles. It reproduces the density, latent heat and phase change of water at a fraction of the cost of all-atom water.
- **Water–metal and water–tool coupling:** Lennard-Jones potentials
- **Backing plate:** the bottom 3 Å of the plates is fixed (rigid backing plate, or anvil). The next 3 Å is held at 300 K (heat sink into the backing plate).
- **Bulk water:** the top 8 Å of water is held at 300 K, standing in for the bulk of the water tank (heat sink into the water)
- **Periodic window:** the box is periodic in x and y, so the plates are infinite in plane. The simulation is a periodic window, about 24 × 12 nm, in the middle of a large submerged plate.
- **Top wall:** a soft wall at the top of the box keeps the water from escaping (a free water surface below a vacuum gap)

---

## Process Parameters

| Quantity                          | Value                  | Equivalent                 |
| --------------------------------- | ---------------------- | -------------------------- |
| Tool rotation `w`                 | 0.05 rad/ps            | shoulder rim speed 125 m/s |
| Plunge speed `vp`                 | 0.10 Å/ps              | 10 m/s                     |
| Plunge depth `dpl`                | 17.5 Å                 | shoulder 0.5 Å into the surface |
| Dwell time `tdw`                  | 20 ps                  |                            |
| Traverse speed `vt`               | 0.25 Å/ps              | 25 m/s                     |
| Traverse length `Ltr`             | 60 Å                   |                            |
| Timestep                          | 1 fs                   |                            |
| Initial temperature               | 300 K                  |                            |

All quantities are in LAMMPS `metal` units (Å, ps, eV). The speeds are orders of magnitude higher than in real welding, which is normal for MD; compare runs with each other (wet vs dry, Al/Al vs Al/Cu) and report trends rather than absolute values.

---

## Model Parameters

| Parameter                         | Value              | Role                                              |
| --------------------------------- | ------------------ | ------------------------------------------------- |
| Box                               | 238.95 × 121.5 × 81 Å | periodic in x, y; free in z                    |
| Plate thickness `zAl`             | 26 Å               | pin reaches about 60 % of the thickness           |
| Al lattice constant               | 4.05 Å             | box length = 59 × 30 whole fcc cells              |
| Cu lattice constant (`cu 1`)      | 3.615 Å            | strained < 0.2 % to fit the box                   |
| Seam position / clearance `sgap`  | y = 60.75 Å / 2 Å  | plate A for y > seam, plate B for y < seam        |
| Pin radius / length               | 10 Å / 15 Å        |                                                   |
| Shoulder radius / height          | 25 Å / 10 Å        |                                                   |
| Tool start position               | x = 70 Å, 2 Å above the surface |                                      |
| Tool offset `yoff`                | 0 Å                | + shifts the pin towards plate A                  |
| Water layer                       | z = 28.5 to 66 Å   | about 13 Å above the shoulder                     |
| Fixed / 300 K plate layers        | z < 3 Å / 3–6 Å    | backing plate                                     |
| 300 K water layer                 | z > 58 Å           | bulk-water reservoir                              |

### Interatomic Potentials

| Pair                    | Style        | Parameters                          | Source / note                          |
| ----------------------- | ------------ | ----------------------------------- | -------------------------------------- |
| Al–Al, Al–Cu, Cu–Cu     | `eam/alloy`  | `AlCu.eam.alloy`                    | LAMMPS `potentials/` folder            |
| water–water (mW)        | `sw`         | `mW.sw`                             | Molinero & Moore, *J. Phys. Chem. B* 113, 4008 (2009) |
| Al–C (tool)             | `morse`      | D0 = 0.28 eV, α = 2.78 Å⁻¹, r0 = 2.20 Å | used in Al/diamond cutting MD     |
| Cu–C (tool)             | `morse`      | D0 = 0.087 eV, α = 5.14 Å⁻¹, r0 = 2.05 Å | used in Cu/diamond cutting MD    |
| Al–water                | `lj/cut`     | ε = 0.0135 eV, σ = 2.90 Å           | starting value, **calibrate**          |
| Cu–water                | `lj/cut`     | ε = 0.0150 eV, σ = 2.75 Å           | starting value, **calibrate**          |
| C–water                 | `lj/cut`     | ε = 0.0050 eV, σ = 3.20 Å           | starting value, **calibrate**          |
| C–C                     | none         |                                     | tool is rigid                          |

The metal–water and tool–water parameters decide how much heat the water removes. Calibrate them against the contact angle of water on the metal or the metal–water interfacial thermal conductance before drawing quantitative conclusions.

---

## Geometry

```
Side view through the tool axis (x-z plane, A):

 z = 76  ------------------------------------------------------------- top wall
          ~~~~~~~~~~~~~~~ 300 K water reservoir (z > 58) ~~~~~~~~~~~~~~~
 z = 66   ~~~~~~~~~~~~~~~~~~~~~~~~~~~ water surface ~~~~~~~~~~~~~~~~~~~~
          ~~~~~~  ______________  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ~~~~~~ |   shoulder   | ~~~~~~   mW water   ~~~~~~~~~~~~~~~~~
          ~~~~~~ |___R = 25_____| ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ~~~~~~~~~~ |pin |  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ~~~~~~~~~~ |R=10|  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 z = 26   ==========================================================  plate surface
          |  plate A (y > seam)  /  plate B (y < seam)               |
          |                                                          |
 z = 6    |- - - - - - - - - - - 300 K layer - - - - - - - - - - - - |
 z = 3    |######################  fixed layer  #####################|  backing plate
 z = 0    ==========================================================
                x = 70 (start)  ------------- traverse 60 A ---------->

Top view (x-y plane):

   y = 121.5 +---------------------------------------------------------+
             |                 plate A: Al (retreating side)           |
   seam      | - - - - - - ( tool )  - - - - - - - - - - - - - - - - - |  y = 60.75
             |                 plate B: Al or Cu (advancing side)      |
   y = 0     +---------------------------------------------------------+
             x = 0                                              x = 238.95

  the tool rotates counter-clockwise (seen from above) and moves in +x:
  the side y < seam (plate B) is the advancing side
```

### System Size

| Component            | Al/Al wet  | Al/Cu wet  | Al/Al dry  | Al/Cu dry  |
| -------------------- | ---------- | ---------- | ---------- | ---------- |
| Plate A, Al (type 1) | 22,243     | 22,243     | 22,243     | 22,243     |
| Plate B (type 4)     | 23,010 Al  | 33,660 Cu  | 23,010 Al  | 33,660 Cu  |
| Tool, C (type 2)     | 4,243      | 4,243      | 4,243      | 4,243      |
| Water, mW (type 3)   | 35,820     | 35,820     | 0          | 0          |
| **Total**            | **85,316** | **95,966** | **49,496** | **60,146** |

---

## Atom Types

| Type | Material                 | Role                                          |
| ---- | ------------------------ | --------------------------------------------- |
| 1    | Al                       | Plate A (retreating side, y > seam)           |
| 2    | C (diamond, rigid)       | Tool: pin + shoulder                          |
| 3    | mW water                 | Water layer                                   |
| 4    | Al (`cu 0`) or Cu (`cu 1`) | Plate B (advancing side, y < seam)          |

---

## Simulation Phases

| Phase         | Steps     | Time    | Description                                                     |
| ------------- | --------- | ------- | --------------------------------------------------------------- |
| Minimisation  | ≤ 2,000   | –       | Removes overlaps from the lattice set-up; tool held fixed       |
| Equilibration | 20,000    | 20 ps   | All mobile atoms thermostatted at 300 K; tool held fixed        |
| Plunge        | 175,000   | 175 ps  | Tool rotates and moves down 17.5 Å at 0.1 Å/ps                  |
| Dwell         | 20,000    | 20 ps   | Tool rotates in place; heat builds up                           |
| Traverse      | 240,000   | 240 ps  | Tool rotates and moves 60 Å along the seam at 0.25 Å/ps         |
| **Total**     | **455,000** | **455 ps** | Timestep 1 fs                                            |

From the start of the plunge on, only the boundary layers are thermostatted (the 300 K plate layer and the 300 K water reservoir); the weld zone heats up freely.

A full Al/Al wet run took **3 h 02 min on 14 MPI ranks** (41.4 timesteps per second, 3.58 ns/day).

---

## Governing Equations

### Tool Kinematics

```
t        = time since the start of the plunge
theta(t) = w t                                              rotation angle
z_shift  = -vp t            for t < tp,   -dpl  after        tp = dpl / vp = 175 ps
x_shift  = vt (t - t1)      for t > t1,    0    before       t1 = tp + tdw = 195 ps

for a tool atom that started at (x0, y0) relative to the axis (xc0, yt):
  rx = x0 - xc0,  ry = y0 - yt
  x  = xc0 + x_shift + rx cos(theta) - ry sin(theta)
  y  = yt            + rx sin(theta) + ry cos(theta)
  z  = z0 + z_shift
  velocities = time derivatives of the above (passed to fix move)
```

### Measured Quantities

```
F_tool   = sum of forces on all tool atoms            (Fz = plunge force, Fx = traverse force)
torque   = sum over tool atoms of (x - xc) fy - (y - yt) fx
mixing   = number of plate-B atoms found within Ly/4 on the plate-A side of the seam

local temperature in a bin (flow removed):
  v_com    = mass-weighted mean velocity of the bin
  T        = sum m |v - v_com|^2 / (3 (N - 1) kB)

unit conversions:  1 eV/A = 1.602 nN,   1 eV = 0.1602 nN nm
```

---

## Group Definitions

```
plA      <- type 1                                 (plate A)
plB      <- type 4                                 (plate B)
metal    <- types 1 and 4
tool     <- type 2                                 (moved by fix move)
water    <- type 3
fixed    <- atoms with z < 3                       (setforce 0, backing plate)
bath     <- metal atoms with 3 < z < 6             (Langevin 300 K)
wtop     <- water atoms with z > zW - 8            (Langevin 300 K, bulk water)
mobile   <- all - fixed - tool                     (integrated with fix nve)
```

---

## Output Files

| File               | Frequency          | Content                                                        |
| ------------------ | ------------------ | -------------------------------------------------------------- |
| `fsw.log`          | every 1,000 steps  | thermo output (see below)                                      |
| `dump.fsw.*`       | every 2,000 steps  | one file per frame: id, type, x, y, z, vx, vy, vz, c_cna, v_Tat (218 frames) |
| `Tmap.dat`         | every 1,000 steps  | surface temperature per 5 Å bin (flow removed)                 |
| `fsw_final.data`   | end of run         | final configuration                                            |

### Thermo Output (welding phase)

| Column     | Quantity                                       | Unit  |
| ---------- | ---------------------------------------------- | ----- |
| `Step`     | Timestep                                       |       |
| `v_t`      | Time since the start of the plunge             | ps    |
| `c_tM`     | Temperature of all metal atoms                 | K     |
| `c_Ft[1]`  | Traverse force Fx on the tool                  | eV/Å  |
| `c_Ft[2]`  | Lateral force Fy on the tool                   | eV/Å  |
| `c_Ft[3]`  | Plunge force Fz on the tool                    | eV/Å  |
| `c_Tq`     | Torque on the tool about its axis              | eV    |
| `c_mix`    | Plate-B atoms carried across the seam          | atoms |
| `v_xsh`    | Tool displacement along the seam               | Å     |
| `v_zsh`    | Tool displacement in z (plunge)                | Å     |

### Per-Atom Columns of `dump.fsw.*`

| Column    | Quantity                                                                  |
| --------- | ------------------------------------------------------------------------- |
| `c_cna`   | Common neighbour analysis: 1 fcc, 2 hcp (stacking faults), 3 bcc, 4 ico, 5 other |
| `v_Tat`   | Instantaneous kinetic temperature of the atom (noisy; average for maps)   |

---

## Repository Structure

```
ufsw_lammps/
|
|-- fsw_underwater.in       # LAMMPS input: plunge, dwell and traverse of the tool (wet or dry)
|-- mW.sw                   # mW water potential (Stillinger-Weber form, metal units)
|-- AlCu.eam.alloy          # Al-Cu EAM potential (copy from the LAMMPS potentials/ folder)
|-- README.md               # this file
|
|-- optional helper files/
|   |-- analyze_fsw.py          # time series, temperature maps, peak temperature, CNA, flow fields
|   |-- render_fsw.py           # renders images or a movie of the dumps without the GUI
|   |-- make_ovito_session.py   # writes ready-made OVITO sessions (atoms, solid surfaces)
|
|-- generated on run/
    |-- fsw.log
    |-- dump.fsw.000000 ... dump.fsw.434000
    |-- Tmap.dat
    |-- fsw_final.data
    |-- fsw_timeseries.png, fsw_Tmaps.png, fsw_peakT.png, fsw_flow.png, fsw_summary.csv
    |-- fsw_atoms.ovito
    |-- fsw_surface.ovito
```

---

## How to Run

### Requirements

- LAMMPS with the **MANYBODY** package (`eam/alloy`, `sw`): <https://www.lammps.org> (tested with the 22 Jul 2025 release)
- An MPI runtime that matches the one LAMMPS was built with, for parallel runs
- Python with numpy, pandas and matplotlib
- OVITO desktop and the `ovito` Python module (**same version** as the desktop app) for visualisation: <https://www.ovito.org>

Check that LAMMPS has the needed styles:

```
lmp -h | grep -E "eam/alloy|sw "
```

### Step 1: Prepare the folder

Put `fsw_underwater.in`, `mW.sw` and `AlCu.eam.alloy` in one folder, for example `D:\extras` (`/mnt/d/extras` in WSL).

### Step 2: Run the weld

With MPI on Ubuntu or WSL:

```
cd /mnt/d/extras
mpirun -np 14 lmp -in fsw_underwater.in -log fsw.log
```

| Option              | Effect                                                   |
| ------------------- | -------------------------------------------------------- |
| `-var wet 0`        | Dry reference run (no water)                             |
| `-var cu 1`         | Dissimilar weld: plate B is Cu                           |
| `-var yoff 3`       | Pin offset 3 Å towards plate A (Al), common for Al/Cu    |
| `-var zAl 30`       | Plate thickness 30 Å (do not go below about 22 Å)        |
| `-var scale 0.5`    | Half-size box for quick tests                            |

Run the wet and the dry case in **separate folders**, so the output files do not overwrite each other:

```
mkdir -p wet dry
cp fsw_underwater.in mW.sw AlCu.eam.alloy wet/ && cp fsw_underwater.in mW.sw AlCu.eam.alloy dry/
cd wet && mpirun -np 14 lmp -in fsw_underwater.in -log fsw.log && cd ..
cd dry && mpirun -np 14 lmp -in fsw_underwater.in -var wet 0 -log fsw.log && cd ..
```

The run is finished when the log ends with `Total wall time`.

### Step 3: Analyse

In the folder of a run:

```
python3 analyze_fsw.py               # add --cu for an Al/Cu run
```

To compare runs (one folder each):

```
python3 analyze_fsw.py wet=/mnt/d/extras/wet dry=/mnt/d/extras/dry
```

| Output                 | Content                                                                     |
| ---------------------- | --------------------------------------------------------------------------- |
| `fsw_timeseries.png`   | Metal temperature, plunge force, traverse force, torque and mixing vs time  |
| `fsw_Tmaps.png`        | Surface temperature maps at the end of the plunge, mid-traverse and the end |
| `fsw_peakT.png`        | Peak surface temperature and fcc/hcp/other fractions of the stirred zone vs time |
| `fsw_flow.png`         | Stirring flow of water and metal around the tool (dwell and traverse)       |
| `fsw_summary.csv`      | Peak temperature, mean traverse forces and torque, final mixing, final fcc fraction |

The temperature maps and the flow fields are built from the dump files: velocities are averaged in small bins and the flow velocity is removed for the temperature, so the random thermal motion does not show up as heat or flow.

### Step 4: Visualise in OVITO

Run in the folder with the dump files. Use **Windows Python** if you use the Windows OVITO app, so the file paths stored in the session are valid:

```
cd /d D:\extras
pip install ovito==<your OVITO version>
python make_ovito_session.py         # add --cu for an Al/Cu run
```

Open a session with **File > Load Session State** (not Load File). The session opens on the last frame.

| Session file          | What it shows                                                                                 |
| --------------------- | --------------------------------------------------------------------------------------------- |
| `fsw_atoms.ovito`     | Atoms as spheres: plate A grey, plate B orange, tool dark, water blue; the front half of the water is removed so the plates and the tool are visible behind it |
| `fsw_surface.ovito`   | Plates, tool and water as solid closed surfaces (water semi-transparent)                     |

Both sessions are set up for the Tachyon renderer with ambient occlusion at 1920 × 1080. Use **Render active viewport** for an image, or set the Render tab to **Complete animation** with an output file such as `fsw.mp4` for a movie.

To render without the GUI:

```
python3 render_fsw.py                          # last frame -> fsw.png
python3 render_fsw.py --water hide             # plates and tool only
python3 render_fsw.py --color temp --tmax 900  # metal coloured by temperature
python3 render_fsw.py --color speed            # metal coloured by speed
python3 render_fsw.py --movie                  # all frames -> fsw.mp4
```

---

## Important: Things That Strongly Affect the Results

| Setting                                      | Effect                                                                    |
| -------------------------------------------- | ------------------------------------------------------------------------- |
| Metal–water LJ ε                             | Sets the heat flow into the water; too weak and the water hardly cools the weld |
| Thermostat only on the boundary layers       | A thermostat on the weld zone would remove the frictional heat artificially |
| Plate thickness vs pin length                | If the pin reaches the fixed or 300 K layers, temperatures and forces are distorted |
| Rotation and traverse speed                  | MD speeds are very high; check how sensitive the results are to `w` and `vt` |
| Water depth                                  | The 300 K water reservoir stands in for the bulk water; a thinner water layer cools more strongly |

---

## What to Look for in Results

### Plunge

The plunge force Fz rises as the pin enters the plates and peaks when the shoulder touches the surface. The metal temperature rises with plastic work under the pin.

### Dwell

The tool rotates in place. Heat builds up under the shoulder; in the wet run part of it flows into the water. `fsw_flow.png` shows the metal swirling around the pin.

### Traverse

The traverse force Fx and the torque settle to a steady level. Material is carried from the advancing side (plate B) around the pin to the retreating side; `c_mix` grows. Behind the tool, the stirred zone cools and partly recrystallises (fcc fraction in `fsw_peakT.png`).

### Wet vs Dry

The main result is the comparison: lower peak temperature and faster cooling under water, a narrower hot zone in `fsw_Tmaps.png`, and the effect on forces, torque, mixing and the defect content of the weld.

### Mixing

Colour by type in OVITO: grey (plate A) and orange (plate B) atoms intermixed in the stirred zone show the material flow across the seam.

---

## Common Errors and Fixes

| Error or symptom                                         | Cause                                                       | Fix                                                                      |
| -------------------------------------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------ |
| `cd: D:\extras: No such file` in Ubuntu                  | Windows path typed in a Linux (WSL) terminal                | Use `cd /mnt/d/extras`                                                   |
| `lmp: command not found`                                 | Different executable name, or LAMMPS not on the PATH        | Try `lmp_mpi`, or add the LAMMPS `bin` folder to the PATH                |
| `Cannot open EAM potential file AlCu.eam.alloy`          | Potential file not in the run folder                        | Copy it from the LAMMPS `potentials/` folder                             |
| `Unrecognized pair style 'sw'` or `'eam/alloy'`          | LAMMPS built without the MANYBODY package                   | Install a LAMMPS that includes MANYBODY (check with `lmp -h`)            |
| mpirun refuses to start 14 processes                     | Fewer cores available                                       | Use fewer ranks, or add `--oversubscribe`                                |
| `Tmap.dat` holds only 6 numbers per output               | Run made with an early version of the input                 | Use `analyze_fsw.py`, which builds the temperature maps from the dumps   |
| Only points in OVITO                                     | Dump files opened directly; default radius 0.5 Å            | Open a session with **File > Load Session State**                        |
| Everything looks blue in a Tachyon render                | Thick water layer hides the plates                          | Use `fsw_atoms.ovito` (front half of the water removed)                  |
| Session cannot find the dump files                       | Session built in WSL (`/mnt/d/...`) and opened in Windows   | Build it with Windows Python, or use **Pick new file** on the file source |
| Session does not open, or version error                  | `ovito` Python module and OVITO app have different versions | `pip install ovito==<exact app version>`                                 |
| No mixing visible                                        | Viewing frame 0, before the tool enters the plates          | Go to the last frame (the sessions open there)                           |
| Qt, `libOpenGL` or Vulkan error when rendering           | No graphics libraries on WSL or a server                    | `sudo apt install libopengl0 libegl1 mesa-vulkan-drivers`                |
| `ModuleNotFoundError: pandas` or `matplotlib`            | Python without the analysis libraries                       | `pip install numpy pandas matplotlib`                                    |

---

## Extending the Model

| Extension                       | What to change                                                                       |
| ------------------------------- | ------------------------------------------------------------------------------------ |
| Dissimilar Al/Cu weld           | `-var cu 1`, optionally `-var yoff 3` to offset the pin into the Al                  |
| Different plate thickness       | `-var zAl <thickness>`; the water and boundary layers follow automatically           |
| Other process speeds            | `w`, `vp`, `vt` in the welding section of the input                                  |
| Longer weld                     | `Ltr` (traverse length); keep the tool away from the box edge in x                   |
| Different tool shape            | `Rp`, `Lp`, `Rs`, `Hs`, or a threaded/tapered pin by changing the tool regions       |
| Other alloys                    | Replace `AlCu.eam.alloy` with an EAM file for your elements and update `pair_coeff`  |
| All-atom water                  | TIP4P/2005 with `kspace_style pppm/tip4p` and `fix shake` (about 20 × the cost)      |
| Denser output for flow analysis | Dump every 500 steps instead of 2,000 (more disk space)                              |

---

## Citation

If you use this code in your research, please cite:

```
@software{mishra_2026_ufsw_lammps,
  author    = {Mishra, Akshansh},
  title     = {Underwater Friction Stir Welding (UFSW) of Aluminium Plates
               Using LAMMPS Molecular Dynamics},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24296019},
  url       = {https://doi.org/10.5281/zenodo.24296019}
}
```

Plain text citation:

> Mishra, A. (2026). *Underwater Friction Stir Welding (UFSW) of Aluminium Plates Using LAMMPS Molecular Dynamics* [Computer software]. Zenodo. <https://doi.org/10.5281/zenodo.24296019>

Please also cite the potentials used: the mW water model (Molinero & Moore, *J. Phys. Chem. B* 113, 4008, 2009) and the source of the EAM file.

---

## Author

**Akshansh Mishra**
GitHub: <https://github.com/akshansh11>

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />
This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

You are free to:

- **Share**: copy and redistribute the material in any medium or format
- **Adapt**: remix, transform, and build upon the material

Under the following terms:

- **Attribution**: You must give appropriate credit to Akshansh Mishra and provide a link to this repository
- **NonCommercial**: You may not use the material for commercial purposes

Copyright 2026 Akshansh Mishra. All rights reserved for commercial use.
