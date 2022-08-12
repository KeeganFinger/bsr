set(BSR_PHOT_COMMON_SOURCES
  zaface.f90		crees.f
)

add_executable(bsr_phot
  ${BSR_PHOT_COMMON_SOURCES}
  bsr_phot.f90
  mod_phot.f90
  sort_photo.f90
  read_arg.f90
  ak_coef.f90
  phot_sec.f90
  read_bound.f90
  read_Ddat.f90
  read_Hdat.f90
  read_Wdat.f90
)
target_link_libraries(bsr_phot zcom ${BLAS_LAPACK_LIBRARIES})

if(MPI_Fortran_FOUND)
  add_executable(bsr_phot_mpi
    ${BSR_PHOT_COMMON_SOURCES}
    bsr_phot_mpi.f90
    mod_phot_mpi.f90
    sort_photo_mpi.f90
    sym_mat1.f90
    read_arg.f90
    br_arg.f90
    ak_coef_mpi.f90
    phot_sec_mpi.f90
    read_bound_mpi.f90
    read_Ddat_mpi.f90
    read_Hdat_mpi.f90
    read_Wdat_mpi.f90
  )
  target_link_libraries(bsr_phot_mpi zcom ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
endif()
