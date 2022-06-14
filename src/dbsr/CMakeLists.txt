add_subdirectory(DBSR_BREIT3) #DONE
add_subdirectory(DBSR_CI3) #DONE
add_subdirectory(DBSR_CONF3) #DONE
add_subdirectory(DBSR_DMAT3) #DONE
add_subdirectory(DBSR_HD3) #DONE
add_subdirectory(DBSR_HF) #DONE
add_subdirectory(DBSR_MCHF) #DONE
add_subdirectory(DBSR_MULT3) #DONE
add_subdirectory(DBSR_POL3) #DONE
add_subdirectory(DBSR_PREP3) #DONE
add_subdirectory(DBSR_UTILS) #DONE
add_subdirectory(DBSR_jj2ls) #DONE

if(${MPI_Fortran_FOUND})
  add_subdirectory(DBSR_BREIT3_MPI)
  add_subdirectory(DBSR_HD3_MPI)
  add_subdirectory(DBSR_MAT3_MPI)
endif()
