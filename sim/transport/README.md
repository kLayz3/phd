# FOOT telescope transport

`foot_transport` is separate from the reaction executable: no fragmentation process, decay, or
hadronic reaction physics is registered. It retains electromagnetic option4,
including multiple scattering and ionisation. The electron production range
cut is 10 mm, as in the existing reaction simulation; this suppresses explicit
delta rays without removing continuous energy loss.

## Build

Use the Geant4 environment from your installation, then configure the existing
simulation build directory:

```sh
cmake -S sim -B sim/build -DGeant4_DIR=/path/to/lib/cmake/Geant4 \
    -DROOT_DIR=/path/to/root/cmake
cmake --build sim/build --target foot_transport -j4
```

The installed `geant4.sh` sets the dataset paths. Use `--threads 8` with an MT-enabled Geant4 build; a serial-only build rejects requests for more than one thread. ROOT Core is needed because
`phy::Nucleus` currently includes the project's geometry/monad headers. The
existing MPhysics, Geometry, FromChars, TextUtils, JSONParser and FilesystemUtils
implementations are compiled directly; neither libStructures nor Go4 is loaded.
The reaction executable and its macro live under `reactions/8b->7be+p/`.
The telescope configuration loader is shared in `include/TransportConfig.hh`
and `src/TransportConfig.cc`.

## Inputs and invocation

```sh
# Check the parsed inputs without constructing Geant4 or writing output.
./sim/build/transport/foot_transport --ion 9C --energy-per-u 700 --check-config

# Example energy only: supply the appropriate kinetic energy at telescope entry.
./sim/build/transport/foot_transport --ion 9C --energy-per-u 700 \
    --events 100000 --threads 1 --seed 12345 --output transport_9C.root
```

`--setup` and `--cov` default to this checkout's `params/foot_setup.json` and
`params/foot_cov.json`. They use the existing preprocessed `ParseJSON` API, once
on the main thread. Isotopes use `phy::Nucleus::get_ion`, with fully stripped
nuclei required. Both the project's tabulated mass and Geant4's particle mass
are recorded; Geant4 retains its own particle definition and transport mass.

- `DSSSD.composition` is repeated for **each of the eight detectors**. Thickness
  is in metres in the input. Use `Si` or `G4_Si` for the single measurement slab;
  other materials use Geant4 NIST material names, e.g. `G4_Al`.
- Composition order is upstream to downstream. Each detector's configured z is
  the silicon midplane. Other stack components are positioned relative to that
  plane; exactly one silicon slab is required in each repeated stack.
- Detector z is `box.det_pos[N] + detector.dz - box.target.dz`, matching FOOT
  coordinates with the target downstream face at zero. Missing detector `dz`
  defaults to zero. Readout rows match `TFOOTHitProc::SetPairConversionMatrices`:
  X/u uses `(cos(angle), sin(angle))`, Y/v uses `(-sin(angle), cos(angle))`.
  The bare measurement is `a*x + b*y - delta_p`. Orientation signs/reflections
  must not be applied again to this already geometric measurement.
- Pair numbers in `foot_cov.json` are one-based. Residual widths are input in
  micrometres, converted to mm, and saved as reference data. Missing isotope
  measurements in a pair are marked unavailable. They are **not** used to smear
  the Geant4 crossings: detector noise is the unknown to be inferred later.

The current transport geometry is a vacuum world and eight axis-aligned square
material stacks, with default half-width 48 mm. Their centres are offset along
the readout axis by `delta_p`; readout rotations affect the measurement, not a
physical tilt. `--half-width-mm` controls this approximate acceptance. Detailed
mechanical outlines, supports, air and target material are not inferred from
unrelated setup keys. Only materials explicitly listed in composition are
included. Exact physical transverse acceptance is not yet a calibration model.

The source starts 1 mm upstream of the first complete stack. `--x-mm`, `--y-mm`,
`--ax` and `--ay` specify its position and slopes **at that entrance plane**.
Independent Gaussian beam spreads are available through `--beam-sigma-mm` and
`--angle-sigma-mrad`. Defaults are an on-axis pencil beam. Broad samples can miss
planes; missing crossings are retained, not silently rejected.

## Output

The ROOT TTree `crossings` contains one row per generated event, with merging
across Geant4 workers. Row order is not event order in multithreaded runs.

- `event`: Geant4 event ID.
- `beam_x_mm`, `beam_y_mm`, `beam_ax`, `beam_ay`: generated entrance state.
- `d0_*` through `d7_*`: physical detector index `N`, **not** pair u/v order.
  Each has `present`, `x_mm`, `y_mm`, `z_mm`, `ax`, `ay`, `kinetic_MeV`, and
  `measurement_mm`. A detector's u/v identity and pair are in the metadata.
- `z_mm` is in FOOT coordinates, not the internally shifted Geant4 world frame.
- Missing values have `present=0` and NaN floating fields. Only the first forward
  crossing of the primary particle is recorded; secondary tracks are not scored.

Silicon is split into two touching half-slabs so its midplane is an exact
geometry boundary. A measurement is taken at the primary's upstream-half to
downstream-half transition. This is a point-crossing approximation, not an
energy-weighted strip centroid. Keep the distinction when modelling response.

`OUTPUT.root.json` records the input file paths, resolved planes and material stack,
reference residuals, isotope masses, seed, source settings, material scope and
step limit. Existing output/metadata files are rejected rather than overwritten.
Parent output directories must already exist.

## Validation and next analysis step

`--vacuum` replaces every detector material by vacuum while retaining the
measurement boundaries. With a fixed source direction, all eight crossings
should lie on the generated line to numerical tolerance. It provides a direct
check of geometry, axes, offsets and detector indexing.

The material step limit defaults to 0.015 mm. It uses `G4StepLimiterPhysics`;
vary `--max-step-mm` to establish convergence before interpreting residual tails.
Production cuts and the simplified material budget also need validation.

For resolution calibration, reuse these truth crossings while varying Gaussian
measurement noise, then use the existing unweighted fit and leave-one-pair-out
residual definitions. Compare against the selected reference widths with the same
core-fitting procedure. The simulation does **not yet** implement detector charge
response, CKF branching/pruning, or the analysis charge gates. Consequently a
spatial-only score cut cannot claim to reproduce `GetScore()`, which also includes
charge residuals. Those selection effects require an explicit downstream response
and reconstruction stage before quoting calibrated intrinsic resolutions.

## Quick plots

From the repository root:

```sh
root -l 'sim/transport/draw_transport.C("transport_9C.root")'
# Optional four-page PDF, also usable in batch mode:
root -l -b -q 'sim/transport/draw_transport.C("transport_9C.root", "transport.pdf")'
```

The canvases show crossing multiplicity, scalar measurements and true x/y at
all eight midplanes, accumulated x/y slope changes relative to the generated
entrance direction, and energy loss between the first and last midplanes.
Detector labels use physical indices d0..d7. The angular labels express slope
differences times 1000, approximating projected angles in mrad near the beam axis.

## Leave-one-pair-out residuals with the analysis binning

Build the ROOT-only analysis target and run it on an existing transport file:

```sh
cmake --build sim/build --target foot_transport_analysis -j4
./sim/build/transport/foot_transport_analysis -f transport_9C.root --section c9
```

Use the ROOT environment matching `ROOT_DIR` in the CMake cache. In particular,
do not combine FairSoft ROOT libraries with a different installation's
`ROOTSYS`/`ROOT_INCLUDE_PATH`; source the matching `thisroot.sh` first.

This produces `transport_9C_residuals.root` and `.pdf`. It reads the source file's
compact `.root.json` summary to construct the calibrated scalar design rows,
then refits the other three pairs for each omitted pair. The PDF contains four
u and four v residual histograms, their Gaussian core fits, and the corresponding
measured core widths from the input summary. All eight crossings must be present.

Binning/settings precedence is:

1. Compiled `foot_track_analysis.cc` defaults: u/v `100,-5,5` in mm, fit side ratio 2.
2. `GLOBAL()` in `quick/hit/cfg/foot_track_analysis/config.cfg`.
3. The selected `SECTION(...)`, using the existing section/argument parsing APIs.
4. Explicit residual-analysis CLI overrides.

The actual carbon section name is **c9**, not 9c. Currently it inherits u/v
`180,-0.6,0.6` from GLOBAL, giving a bin width of 6.67 micrometres. Other config
arguments such as files and charge gates are ignored here. `--no-config` uses
compiled defaults; `--section ''` reads only GLOBAL. Gaussian fitting uses the
same `GaussFitMax` routine and four iterations as `foot_track_analysis.cc`.

```sh
# Change only u binning; v still comes from the selected configuration.
./sim/build/transport/foot_transport_analysis -f transport_9C.root \
    --binning-dx=240,-0.6,0.6 -o finer_residuals.root

# Optional trial detector noise: physical detector order d0..d7, in micrometres.
./sim/build/transport/foot_transport_analysis -f transport_9C.root \
    --sigma-um=30,30,30,30,30,30,30,30 -o trial_residuals.root
```

Noise defaults to zero, so the first comparison is scattering-only. Trial noise
is added to measurements before fitting; transport is not rerun. The ROOT file
also stores the residual response R, variance response T, full empirical residual
covariance (including values outside the histogram range), and a small analysis
summary. The mapping satisfies `residual = R * (measurement + offset)` and
`Var(residual)[i] = sum_j T[i,j] * detector_variance[j] + scattering_variance[i]`.
The variance response rank is reported as a diagnostic, not a guarantee that all
eight detector variances can be estimated accurately from core widths.

This is a forward residual comparison, **not yet an automatic resolution fit**.
It applies no CKF pruning, charge response/gates or final track-score cut, so
matching the binning and Gaussian fitting does not yet match real-data selection.
Any detector-resolution inference still needs the selection-bias validation
specified in the uncertainty plan.

A larger transport sample can be requested with `--events 150000`. The analysis
reads every input event; it does not cap the sample at 10k or 150k.
