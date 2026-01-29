#### Analysis for my PhD
All the code presented here is meant for (personal) use. 
Feel free to use it or fork it.

Remember:
`git submodule update --init --recursive`

Ceres library isn't needed for main analysis yet, just as a routine for GaussFitter.
For `ceres`:

Its dependencies are submoduled!
To build it:
```
cmake -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=ON
```

See how its linked properly in `GaussFitter.hxx` file. It's a linking nightmare.

--𝒦𝓁𝒶𝓎𝓏𝓮
