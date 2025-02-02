PROJECT(bsr)
CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
INCLUDE(CheckLanguage)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mods)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

ENABLE_LANGUAGE(Fortran)
CHECK_LANGUAGE(CUDA)

if (CMAKE_CUDA_COMPILER)
  MESSAGE(STATUS "Enabling CUDA")
  ENABLE_LANGUAGE(CUDA)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -stdpar -gpu=managed -cudalib=cusolver,cublas")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  MESSAGE(STATUS "Setting build type to 'Debug' as none was specified.")
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build." FORCE)
endif()

GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  MESSAGE(STATUS "Compiling using gfortran")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-automatic -fcray-pointer -fno-f2c -O2 -fallow-argument-mismatch")
  SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  MESSAGE(STATUS "Compiling using pgfortran")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-Mpreprocess -Munroll -Msave -Mcray=pointer -O2")
  SET(CMAKE_Fortran_FLAGS_DEBUG "-Mpreprocess -O0 -g")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  MESSAGE(STATUS "Compiling using ifort")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-fpp -unroll-aggressive -save -safe-cray-ptr -O2 -falias")
  SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -mkl")
endif() 

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

FIND_PACKAGE(MPI)
if (${MPI_Fortran_FOUND})
  MESSAGE(STATUS "MPI found - building parallel executables")
endif()

SET(FC_MALLOC LINUX)
SET(BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${SCALAPACK_LIBRARIES})

ADD_SUBDIRECTORY(libs)
ADD_SUBDIRECTORY(bsr)
ADD_SUBDIRECTORY(dbsr)
