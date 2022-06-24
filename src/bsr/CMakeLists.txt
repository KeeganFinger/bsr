add_subdirectory(BSR_BMAT) #DONE
add_subdirectory(BSR_BREIT) #DONE
add_subdirectory(BSR_CONF) #DONE
#add_subdirectory(BSR_DMAT3)
add_subdirectory(BSR_DMAT4) #DONE
#add_subdirectory(BSR_HD3-check)
add_subdirectory(BSR_HD4) #DONE
#add_subdirectory(BSR_MAT3)
#add_subdirectory(BSR_MAT3_MPI)
add_subdirectory(BSR_MAT4) #DONE
add_subdirectory(BSR_MULT) #DONE
add_subdirectory(BSR_PHOT3) #DONE
if(MPI_Fortran_FOUND)
add_subdirectory(BSR_PHOT3_MPI) 
endif()
add_subdirectory(BSR_POL3) #DONE
add_subdirectory(BSR_PREP3) #DONE
add_subdirectory(BSR_RECOUP4) #DONE
add_subdirectory(BSR_UTILS) #DONE
if(MPI_Fortran_FOUND)
  add_subdirectory(BSR_HDB4) #SO MANY ERRORS
endif()
