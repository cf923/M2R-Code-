#!/bin/sh
# run in the same working directory as the C++ file.
# installs dependencies, compiles files, runs interfacing files. Do not run this each time you want to run the code, as this may install a new copy of eigen every time.
# this currently has absolutely no safeguarding. probably breaks.

git clone https://gitlab.com/libeigen/eigen.git
g++ -I eigen/Eigen TISE_solver.cpp -o TISE_solver
python3 interface_TISE_solver.py
echo "if a window popped up with some working plots, then verything resolved, compiled, and ran correctly."