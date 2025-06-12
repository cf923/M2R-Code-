This directory contains numeric solvers for the Classical PT-symmetry and Schrödingers equation.


### Dependencies:

The Eigen C++ library must be downloaded, and usable either in your compiler's include path, or just in the same location as your file. Then compile **including the Eigen folder in the eigen library**. The official Eigen site is, at the time of writing this, down, so clone the GitLab repo:

<https://gitlab.com/libeigen/eigen.git>

For example, to compile the TISE_solver.cpp file to a binary:

all commands ran in the same working directory where the TISE_solver.cpp file is located:

### Unix-based:
``` Bash
g++ -I absolute/path/to/Eigen/directory TISE_solver.cpp -o TISE_solver 
```

### Windows:
using WSL/MinGW-w64/Cygwin:
``` Bash
g++ -I absolute\path\to\Eigen\folder TISE_solver.cpp -o TISE_solver.exe 
```
(with g++ in the PATH)

Alternatively, use VScode.
### macOS:
``` Bash
clang++ -I absolute_path_to_Eigen_directory TISE_solver.cpp -o TISE_solver
```
(g++ will probably be symlinked to clang++, so the Unix command will probably work exactly)

Alternatively, use Xcode.

A Python installation with numpy and matplotlib must exist for all interfacing and plots. Modify the code to have the correct files in any subprocess.run() lines.

