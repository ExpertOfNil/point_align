# Point Alignment

Currently just a library providing a point alignment solver that uses Berthold K.P. Horn's
["Closed-Form Aolution of Absolute Orientation Using Unit Quaternions"](Absolute_Orientation.pdf).
A static binary of OpenBLAS is provided for facilitating eigenvalue and eigenvector solving.
However, an implementation of Ferrari's method is included.

## Setup

1. Clone the repository
1. `cd point_align`
1. `cmake -B build` (add `-DUSE_OPENBLAS=OFF` to disable OpenBLAS and use builtin Ferrari method)
    * Note: the default build mode is `Debug`.  Use `-DCMAKE_BUILD_TYPE=Release` or
    `-DCMAKE_BUILD_TYPE=RelWithDebInfo` for optimizations.
1. Static and dynamic library binaries are compiled to `lib`
1. Examples are compiled to `bin/examples`
1. Tests are compiled to `bin/tests`
