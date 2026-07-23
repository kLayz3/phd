Here I stash some intemediate calculations, that do not convert between ROOT files,
but have become too clunky to fit into a simple ROOT macro.

Structure is the following:
The base program `run` parses thru a config file and delegates the corresponding cmd-line args inside of 
individual named `SECTION(..)` blocks to whichever underlying program in `cal/` or `hit/` or `map/`.

Config files are first parsed thru a GCC preproc, so comment notations are fine.
