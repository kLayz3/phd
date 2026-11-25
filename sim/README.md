# FOOT simulations

Source `sim/sourceme.sh` in the build shell, then build from the repository root.

```sh
cmake -S sim -B sim/build
cmake --build sim/build -j4
sim/build/reactions/reactions --which=1,2 --probability 0.3,0.7 --events 10000 --save=root,pdf --info=boron
```

## Prompt decays

`--list` defaults to [reaction_list.json](reactions/reaction_list.json). The secondary
beam starts just upstream of Be at the listed kinetic energy in MeV/u. The
target file path is relative to the reaction list. Use `--energy-per-u` or
`--setup` for an explicit override. The Be downstream face defines FOOT z = 0.

Each incident particle is assigned at most one channel. Probabilities refer to
incident particles and any remainder gives a beam particle that traverses the setup
without the modeled decay. For a selected channel the production depth is uniform
along z. Geant4 transports the beam to that depth. The excited mother inherits
the beam direction and **total kinetic energy** there. Excitations are sampled in
keV from a Breit–Wigner within five widths of the nominal value and above the
daughter mass threshold. A zero width selects the nominal energy. Decay is
prompt, isotropic two-body or direct unweighted many-body phase space.
The two-proton model omits intermediate resonances and final-state interactions.

The serial simulation transports daughters in Be and in the same material stacks
as `foot_transport`. Each silicon midplane supplies one unsmeared strip-axis
projection. A daughter needs all eight projections for its Eigen straight-line fit.
The event vertex and proton-emission ρ use the same geometry and mass-weighted
angular routines as `quick/hit/physics.cc`. Tracks are associated by simulated
daughter identity. Scintillator gates, charge response and combinatorial track
scores require detector-response models.

Use `--save=root` to write the FOOT canvases individually under
`autosave/simulated/reactions/root`. Use `--info` for additional output directory
components.

## Scattering-only transport

See [the transport notes](transport/README.md) for `foot_transport` and
`analyse_transport`, geometry inputs, plots and resolution studies.
