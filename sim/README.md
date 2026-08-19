# Disclaimer: AI written. I have no clue what it's doing.

## 9C + 9Be toy fragmentation transport in Geant4

Minimal Geant4 application for the controlled model

    9C --(EM transport in Be)--> stochastic reaction vertex
       -> 8B* (effective production model)
       -> p + 7Be
       --(EM transport in remaining Be)--> target exit

The purpose is to isolate the effect of reaction depth, dE/dx and multiple
Coulomb scattering on the p--7Be opening angle.

## What is physical vs deliberately simplified

Geant4 handles:
- 9C, proton and 7Be transport;
- ionisation / energy loss;
- energy-loss fluctuations;
- multiple Coulomb scattering;
- geometry and the stochastic competition between transport steps and the
  custom discrete reaction process.

The custom model currently assumes:
- constant TOTAL 9C reaction/removal cross section = 1 barn (placeholder);
- every nuclear reaction that occurs is forced into the selected 8B* channel;
- 8B* keeps the instantaneous 9C velocity at the reaction point;
- 50/50 population of the 0.770 and 2.320 MeV states;
- truncated Breit-Wigner/Cauchy line shapes with widths 35.6 and 350 keV;
- isotropic 8B* -> p + 7Be decay;
- the unobserved fragmentation system X is not transported;
- no secondary hadronic reactions after the selected fragmentation.

These assumptions are intentionally centralized in `include/Config.hh` and
`src/C9FragmentationProcess.cc`.

## Build

```bash
mkdir build
cd build
cmake ..
cmake --build . -j
./c9be9 run.mac
```

If Geant4 is installed under a non-standard prefix, source its setup script or
pass `-DGeant4_DIR=...` to CMake.

The run writes `c9frag.root` using Geant4's analysis manager.

## Ntuple columns

`events` contains:
- `event`
- `reacted`
- `resonance` (1 = 0.77 MeV, 2 = 2.32 MeV)
- `reaction_z_mm`
- `depth_mg_cm2`
- `c9_T_per_u_MeV`
- `Ex_MeV`
- `theta_vertex_deg`
- `theta_exit_deg`
- `both_exited`
- vertex and exit momentum magnitudes for p and 7Be

The first useful ROOT plots are, for example:

```cpp
ROOT::RDataFrame df("events", "c9frag.root");
auto d = df.Filter("reacted && both_exited");
d.Histo2D({"h", ";reaction depth [mg/cm^{2}];#theta_{exit} [deg]",
           100, 0, 7500, 180, 0, 6},
          "depth_mg_cm2", "theta_exit_deg")->Draw("COLZ");
```

and compare `theta_vertex_deg` against `theta_exit_deg` separately for the two
resonances.

## First validation checks

1. Disable transport conceptually by selecting reactions near the target exit:
   `depth_mg_cm2 > 7400`; theta_exit should approach theta_vertex.
2. Plot reaction depth for reacted events. With a constant cross section it
   should follow a truncated exponential in traversed path length.
3. Split by resonance. The 2.32 MeV state should have the broader/larger
   generator-level p--7Be opening-angle scale.
4. Vary `sigmaTotal`: this should change the reaction-depth distribution, not
   the two-body decay kinematics at a fixed beam energy.
5. Later vary an explicit target step limit and verify the opening-angle result
   is numerically stable before trusting MCS tails.

## Geant4 11.4 custom-process registration

`C9FragmentationProcess` is registered directly with the `G4GenericIon`
`G4ProcessManager` using `AddDiscreteProcess()`.  Do not replace this with
`G4PhysicsListHelper::RegisterProcess()` unless the custom process is also given
a Geant4 process type/subtype that exists in the helper ordering table; otherwise
Geant4 11.4 aborts during initialization with `Run0108: No Matching process
Type/SubType`.


## Delta-electron suppression

The physics list keeps the default production range cut at 0.1 mm but raises
only the secondary-electron (`e-`) production range cut to 10 mm.  This
strongly suppresses explicit delta-ray tracks while leaving sub-threshold
ionisation energy transfer in the continuous stopping power.  Tune
`cfg::deltaElectronRangeCut` in `include/Config.hh`.

This is a **production cut**, not a tracking cut.  It does not kill the decay
proton or the 7Be fragment.

## Multithreading

`main.cc` now uses `G4RunManagerFactory::CreateRunManager()` without forcing a
serial run manager.  When the Geant4 installation has MT support, the
`/run/numberOfThreads 8` line in `run.mac` selects eight workers.  It must occur
before `/run/initialize`.  ROOT ntuple merging is enabled, so worker ntuples
are merged into the main output file.

You can alternatively force the worker count from the shell, for example:

```sh
G4FORCENUMBEROFTHREADS=16 ./c9be9 run.mac
```

The environment variable overrides the run-manager/UI thread count in an MT
build.
