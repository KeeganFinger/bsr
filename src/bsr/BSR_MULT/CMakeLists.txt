#add_subdirectory(mult3)
#add_subdirectory(mult3_tab)
add_subdirectory(mult4)
add_subdirectory(simple2)

if(MPI_Fortran_FOUND)
  add_subdirectory(mult4_mpi)
endif()
