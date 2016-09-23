opt="-Wextra -pedantic -g -O2 "

cmake .. \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_C_FLAGS="$opt" \
-DCMAKE_CXX_FLAGS="$opt" \
-DCMAKE_EXE_LINKER_FLAGS="-ldl $opt" \
-DCMAKE_INSTALL_PREFIX=$PWD/install_nothread/ \
\
-DSCOREC_PREFIX=/path/to/SCOREC/core/install/ \
\
-DPHASTA_INCOMPRESSIBLE=ON \
-DPHASTA_COMPRESSIBLE=ON \
-DPHASTA_USE_SVLS=OFF \
-DPHASTA_USE_PETSC=OFF \
-DLESLIB=/path/to/libles.a \
-DPHASTA_TESTING=ON \
-DCASES=/path/to/phastaChefTests/ \
..

make
make test
