#!/usr/bin/env bash
# Generate clangd's compilation database with the simulation's FairSoft toolchain.
set -e

sim_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$sim_dir/sourceme.sh"

cmake -S "$sim_dir" -B "$sim_dir/.clangd-build" \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
	-DROOT_DIR="$(root-config --prefix)/share/root/cmake" \
	-DGeant4_DIR="$(geant4-config --prefix)/lib/cmake/Geant4"

cat > "$sim_dir/.clangd" <<'YAML'
CompileFlags:
  CompilationDatabase: .clangd-build
YAML
