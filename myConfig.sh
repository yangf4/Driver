#!/bin/bash -e

opt="-Wextra -pedantic -g -O0 "

cmake .. \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_BUILD_TYPE=Debug \
-DCMAKE_C_FLAGS="$opt" \
-DCMAKE_CXX_FLAGS="$opt" \
-DCMAKE_EXE_LINKER_FLAGS="-ldl $opt" \
-DENABLE_THREADS=OFF \
-DSIM_MPI="mpich3.1.2" \
-DSIM_PARASOLID=ON \
-DPCU_COMPRESS=ON \
-DENABLE_ZOLTAN=ON \
-DIS_TESTING=True \
-DMESHES=/fasttmp/yangf4/phastaChef/meshes/ \
-DCORE_SRC_DIR=/fasttmp/yangf4/fork_chef/core-sim \
-DPHASTA_SRC_DIR=/fasttmp/yangf4/git_phasta/phaseChangePhasta \
-DPHASTA_INCOMPRESSIBLE=OFF \
-DPHASTA_COMPRESSIBLE=ON \
-DLESLIB=LIBLES \
-DCASES=/lore/shamse/codes/scorec_phasta/phastaChefTests \
-DPHASTA_TESTING=ON \
..
