#### Analysis for my PhD
All the code presented here is meant for (personal) use. 
Feel free to use it or fork it.
``
git clone --recurse-submodules
``

Only external dependency should be [https://www.boost.org/doc/user-guide/getting-started.html](boost).
Note that few of the subroutines from `util/` won't compile on non-POSIX systems.

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

Python3 only needed for the pretty histogramming API. On GSI cluster is a nightmare to due to various
wrong pip versions.
In this case:
```
python3 -m venv --without-pip .venv
curl -sS https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py
.venv/bin/python /tmp/get-pip.py
.venv/bin/python -m pip install numpy matplotlib
set -gx PYTHONPATH "$PWD/.venv/lib/python3.13/site-packages"
```

Later just `source .venv/bin/activate` if working in bash/zsh, otherwise `activate.fish`

--𝒦𝓁𝒶𝓎𝓏𝓮
