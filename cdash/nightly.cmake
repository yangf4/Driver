cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT ON)
SET(CTEST_TEST_TYPE Nightly)

set(CTEST_SITE             "pachisi.scorec.rpi.edu" )
set(CTEST_DASHBOARD_ROOT   "/lore/cwsmith/cdash" )
set(CTEST_CMAKE_GENERATOR  "Unix Makefiles" )
set(CTEST_BUILD_CONFIGURATION  RelWithDebInfo)

set(CTEST_PROJECT_NAME "PhastaChef")
set(CTEST_SOURCE_NAME repos)
set(CTEST_BUILD_NAME  "linux-gcc-${CTEST_BUILD_CONFIGURATION}")
set(CTEST_BINARY_NAME build)

set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif()
if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
endif()

configure_file(${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake
               ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake COPYONLY)

set(CTEST_NIGHTLY_START_TIME "03:00:00 EST")
set(CTEST_BUILD_FLAGS -j8)

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=phastaChef")
set(CTEST_DROP_SITE_CDASH TRUE)

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

ctest_start(${CTEST_TEST_TYPE})

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit PhastaChef Project.xml!")
  endif()
endif()

macro(submit_part subproject_name part)
  if(CTEST_DO_SUBMIT)
    ctest_submit(PARTS ${part} RETURN_VALUE HAD_ERROR)
    if(HAD_ERROR)
      message(FATAL_ERROR "Cannot submit ${subproject_name} ${part} results!")
    endif()
  endif()
endmacro()

macro(build_subproject subproject_name config_opts)
  set_property(GLOBAL PROPERTY SubProject ${subproject_name})
  set_property(GLOBAL PROPERTY Label ${subproject_name})

  setup_repo(${subproject_name} "git@github.com:PHASTA/phastaChef.git")
  if( ${subproject_name} STREQUAL "phastaChef-sim" )
    setup_repo(core-sim "git@github.com:SCOREC/core-sim.git")
  elseif( ${subproject_name} STREQUAL "phastaChef" )
    setup_repo(core "git@github.com:SCOREC/core.git")
  endif()
  setup_repo(phasta "git@github.com:PHASTA/phasta.git")

  if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}/${subproject_name}")
    file(MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/${subproject_name})
  endif()

  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}" 
    SOURCE "${CTEST_SOURCE_DIRECTORY}/${subproject_name}"
    OPTIONS "${config_opts}"
    RETURN_VALUE HAD_ERROR)
 
  submit_part(${subproject_name} Configure)

  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}"
      RETURN_VALUE HAD_ERROR)

  submit_part(${subproject_name} Build)
endmacro()

macro(test_subproject subproject_name)
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}")

  submit_part(${subproject_name} Test)
endmacro()

macro(setup_repo repo_name repo_url)
  if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/${repo_name}")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
      clone ${repo_url} ${CTEST_SOURCE_DIRECTORY}/${repo_name}
      RESULT_VARIABLE HAD_ERROR)
    if(HAD_ERROR)
      message(FATAL_ERROR "Cannot checkout ${repo_name} repository!")
    endif()
  endif()

  ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/${repo_name}" RETURN_VALUE count)
  message("Found ${count} changed files")
  submit_part(${repo_name} "Update")
endmacro(setup_repo)

set(flags "-O2 -g -Wall -Wextra")
SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_Fortran_COMPILER=gfortran"
  "-DCMAKE_C_FLAGS=${flags}"
  "-DCMAKE_CXX_FLAGS=${flags}"
  "-DCMAKE_Fortran_FLAGS=${flags}"
  "-DCMAKE_EXE_LINKER_FLAGS=${flags}"
  "-DCORE_SRC_DIR=${CTEST_SOURCE_DIRECTORY}/core-sim" ##change this to core for a non-sim build
  "-DSIM_MPI=mpich3.1.2"
  "-DPCU_COMPRESS=ON"
  "-DENABLE_THREADS=OFF"
  "-DENABLE_ZOLTAN=ON"
  "-DIS_TESTING=True"
  "-DMESHES=/lore/cwsmith/develop/streamInMem/meshes"
  "-DPHASTA_SRC_DIR=${CTEST_SOURCE_DIRECTORY}/phasta"
  "-DPHASTA_INCOMPRESSIBLE=ON"
  "-DPHASTA_COMPRESSIBLE=ON"
  "-DPHASTA_TESTING=ON"
  "-DACUSOLVE_LIB=/users/cwsmith/develop/libLes/libles_gcc_nolic.a"
  "-DCASES=/lore/cwsmith/develop/streamInMem/phastaChefTests"
)
build_subproject(phastaChef-sim "${config_opts}")
test_subproject(phastaChef-sim)

message("DONE")