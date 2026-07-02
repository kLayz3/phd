Here I stash some intemediate calculations, that do not convert between ROOT files,
but have become too clunky to fit into a simple ROOT macro.

Structure is the following:
Whatever would represent a ROOT macro, stuff into it's own little function and namespace. Similar to before.
Then the main program will switch based on the first argument into whatever little function it catches.
