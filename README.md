# LSCSimulation

Geant4-based Monte Carlo simulation for the LSC (Liquid Scintillator Counter) Prototype detector.

## Modules

| Module | Description |
|--------|-------------|
| `MCObjs` | ROOT-based MC data objects (tracks, hits, scintillation) |
| `GLG4Sim` | Low-level Geant4 simulation framework (GLG4) |
| `LSCSim` | Detector construction, physics, and ROOT output |
| `TrgSim` | PMT signal and trigger simulation |

## Dependencies

- **Geant4** — set `G4BASE` to the Geant4 installation prefix
- **ROOT** — `root-config` must be on `PATH`
- **CMake** ≥ 3.15

## Build

```bash
git clone https://github.com/cup-software2018/LSCSimulation.git
cd LSCSimulation

cmake -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j$(nproc)
cmake --install build
```

If Geant4 is not found automatically, set `G4BASE` before configuring:

```bash
export G4BASE=/path/to/geant4
cmake -B build -DCMAKE_INSTALL_PREFIX=install
```

## Environment Setup

After installation, **source the setup script before every session**:

```bash
source install/setup_lscsim.sh
```

This script — generated automatically by CMake and installed to the install prefix root — sets:

| Variable | Value |
|----------|-------|
| `PATH` | `install/bin` prepended |
| `LD_LIBRARY_PATH` | `install/lib64` prepended |
| `ROOT_INCLUDE_PATH` | `install/include` prepended |

Without sourcing this script, `lscsim` will not be found and the shared libraries will fail to load. You may want to add it to your shell profile:

```bash
echo "source /path/to/install/setup_lscsim.sh" >> ~/.bashrc
```

## Running

After sourcing the setup script, `lscsim` is available directly:

```
Usage: lscsim [-n events] [-o output.root] [-f macro]
              [-g geometry_data] [-p pmtpos_data] [-m material_data]
              [-v]

  -n  number of events to simulate
  -o  output ROOT file (default: simout.root)
  -f  Geant4 macro file
  -g  geometry data file
  -p  PMT position data file
  -m  material data file
  -v  launch interactive visualization session
```

Data files are in `LSCSim/data/` and example macros in `LSCSim/mac/` of the source tree, and are also copied to `install/share/` on install.

### Example

```bash
source install/setup_lscsim.sh

lscsim -n 1000 -o out.root \
    -f /path/to/LSCSim/mac/default_prototype.mac \
    -g /path/to/LSCSim/data/geometry_prototype.dat \
    -p /path/to/LSCSim/data/pmtpos_prototype.dat \
    -m /path/to/LSCSim/data/materials.dat
```

The macro `default_prototype.mac` uses relative paths (`data/...`), so when using it directly the working directory must contain a `data/` directory or symlink pointing to `LSCSim/data/`.

## Detector

The simulation implements the **Prototype** detector:

- Cylindrical stainless steel buffer tank filled with water
- Inner acrylic target vessel filled with LAB-based liquid scintillator
- 10-inch PMTs positioned from `pmtpos_prototype.dat`
- Optional Teflon reflector cylinder (controlled by `reflector_on` in geometry file)

### Geant4 UI Commands

| Command | Description |
|---------|-------------|
| `/LSC/det/geometrydata <file>` | Geometry parameter file |
| `/LSC/det/materialdata <file>` | Material property file |
| `/LSC/det/pmtposdata <file>` | PMT position file |
| `/LSC/det/geomcheck <0/1>` | Enable geometry overlap check |
| `/LSC/ROOT/savetrackopt <0–4>` | Track saving level (0 = off) |
| `/LSC/ROOT/savestepopt <0/1>` | Save per-step data |
| `/LSC/ROOT/savehitphoton <0/1>` | Save individual photon hits per PMT |
| `/LSC/ROOT/savescintstep <0/1>` | Save scintillation step data |
| `/LSC/Scintillation/scinton <0/1>` | Enable scintillation process |
| `/LSC/Cerenkov/cerenkovon <0/1>` | Enable Cherenkov process |

## Output

The output ROOT file contains a `TTree` with one entry per event:

| Branch | Type | Description |
|--------|------|-------------|
| `PMTData` | `MCPMTData` | PMT hit collection (sorted by PMT ID) |
| `ScintData` | `MCScintData` | Scintillation energy deposition by volume |
| `TrackData` | `MCTrackData` | MC track information |
| `PrimaryData` | `MCPrimaryData` | Primary particle information |
| `EventInfo` | `MCEventInfo` | Run/event metadata |
