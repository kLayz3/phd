Here I stash some intemediate calculations, that do not convert between ROOT files,
but have become too clunky to fit into a simple ROOT macro.

### Configuration
Main program is the `run`. It simply parses out a configuration file and sends it as a sequence of
command-line-arguments inside of individual named `SECTION(...)` to the underlying programs in 
`cal/` or `hit/` or `map/` via exec(1).

Config files are first parsed thru a GCC preproc, so that the comment notations, \#defines and \#includes
are accepted and parsed out.

Successive identical arguments, even for lists, always replace the previously given, never append.

#### Example
```
-a,     --arg INT[,INT...]
```
and program is run as:
```
$PATH/prog -a 1,2,3 --arg=4,5,6
```
delegates `[4,5,6]` to the input vector.

#### Authoritative arguments
Appending the value with `!` character marks the parse value as authoritative, mean that 
subsequent identical arguments, even if themselves authoritative, cannot replace this value.

```
$PATH/prog -a !1,2,3 --arg=4,5,6
```
delegates `[1,2,3]` to the input vector.


