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
-DPCU_COMPRESS=ON \
-DENABLE_ZOLTAN=ON \
-DIS_TESTING=True \
-DMESHES=/fasttmp/yangf4/meshes/ \
\
-DPHASTA_INCOMPRESSIBLE=OFF \
-DPHASTA_COMPRESSIBLE=ON \
-DPHASTA_USE_SVLS=OFF \
-DPHASTA_USE_PETSC=OFF \
-DLESLIB=$LIBLES \
-DPHASTA_TESTING=ON \
-DCASES=/fasttmp/yangf4/phastaChefTests/ \
..

