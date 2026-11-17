# FOOT simulations

```text
sim/
  include/                 shared beam/geometry settings and interfaces
  src/                     shared geometry, beam generator and JSON loader
  transport/               telescope transport and residual-analysis programs
  reactions/
    CMakeLists.txt
    8b->7be+p/             reaction mechanism, parameters, event output and plots
  build/
  CMakeLists.txt
```

`include/Config.hh` contains the incoming beam, target dimensions/density,
world dimensions and EM production cuts used by the reaction simulations.
`include/TransportConfig.hh` and `src/TransportConfig.cc` load the calibrated
FOOT detector geometry and material composition from the existing JSON files.
The transport driver retains its CLI-configured beam and JSON geometry.

Each reaction owns its decay/production model, cross section, resonance widths,
Q-value/threshold calculation, and event/output schema. These belong under its
reaction directory, not in the common geometry configuration.

## Build

After setting up the matching Geant4 and ROOT environments:

```sh
cmake -S sim -B sim/build
cmake --build sim/build -j4
```

Use `-DGeant4_DIR=...` and `-DROOT_DIR=...` on configuration if needed.
Build targets and resulting executables are:

| CMake target | Executable |
|---|---|
| `reaction_8b_to_7be_p` | `sim/build/reactions/8b->7be+p/run` |
| `foot_transport` | `sim/build/transport/foot_transport` |
| `foot_transport_analysis` | `sim/build/transport/foot_transport_analysis` |

Quote reaction paths in the shell: `>` is a redirection operator.
The reaction macro is copied alongside its executable and is also its default
when no argument is supplied. Output ROOT files are written in the current
working directory, so a typical reaction invocation is:

```sh
cd 'sim/build/reactions/8b->7be+p'
./run run.mac
```

The default reaction macro works with serial Geant4. For an MT-enabled build,
uncomment `/run/numberOfThreads 8` before `/run/initialize`.

See [the reaction notes](reactions/8b->7be+p/README.md) for its physics assumptions,
and [the transport notes](transport/README.md) for telescope inputs, truth output,
plots and resolution studies. Existing files in an old build directory may
include obsolete executables; the paths above are the current targets.
