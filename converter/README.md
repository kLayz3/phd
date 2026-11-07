Convert ROOT TTree files using ROOT's RNTupleImporter (ROOT 6.36).

Based on [Jakob Blomer's ATLAS conversion example](https://github.com/jblomer/iotools/blob/master/gen_atlas.cxx)
from iotools, adapted here for more general split Go4/ucesb trees.

`gen_rntuple` imports ordinary supported trees.
ROOT reads the original event using the streamer information stored in the input.

```bash
make -C converter
./converter/gen_rntuple_frs -i input.root -c zstd
./converter/gen_rntuple_frs -i input.root -o /destination -c lz4 -n 1000
./converter/gen_rntuple_frs -i test/testfile.root -t UnpackxTree -c zstd \
    -l ./test/libStructures.so
```

Use `ROOT_CONFIG=/path/to/root-config` with make if ROOT is not on PATH.
Without `-t`, the converter selects the sole TTree in the file, searching subdirectories
as well; use `-t tree/path` to select explicitly.
The default output is
`<input-directory>/rn/rn_<input-filename>`; `-o` selects the output directory
directly. Existing output files are rejected. `-n` optionally limits the entry count.

Base-class branches are deliberately omitted. This is because Go4's `TGo4EventElement`
type down the line inherits from TNamed which is an unhandled case by the default
ROOT converter. Hehe. Every other leaf field is considered:

- Primitive scalar/fixed-array leaves are imported using the stored file schema;
  their enclosing event class does not need a compiled dictionary.
- Object fields and fixed object arrays are imported when a dictionary is available
  and ROOT's RNTuple field factory supports their type. Repeat `-l library.so` to
  load dictionaries before opening the input. Load dependencies before libraries
  that require them.
- Other fields are skipped with their branch name, type and reason. This includes
  unknown helper types accidentally persisted without `//!`; no type-name exclusion
  list is needed. Raw object pointers, variable-length leaf arrays, C-string leaves
  and multi-leaf terminal branches are skipped by the adapter.
  Reason: users should grow up and not use C-api in year 2026.

A missing dictionary for the enclosing type can still produce a ROOT warning;
this does not prevent importing its supported fields. Read failures remain fatal.

ROOT's importer replaces dots in branch names with underscores, for example
`FRSSortEvent.EventFlag` becomes `FRSSortEvent_EventFlag`. Fixed arrays retain
all elements in their original order; ROOT 6.36 imports multidimensional C arrays
as flat `std::array<T, N>` fields, where N is the total number of elements.
Fixed arrays of dictionary-backed objects are also retained as `std::array`.
Histograms and other objects outside the input tree are not copied.
