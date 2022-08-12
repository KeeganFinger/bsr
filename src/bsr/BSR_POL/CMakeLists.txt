add_executable(bsr_pol
  mod_pol.f90              bsr_pol.f90              read_bsw.f90
  read_data.f90            read_arg.f90             read_bsrmat.f90
  read_dipmat.f90          solv_mat.f90             read_nortb.f90
  inf_bsr_pol.f90)#	   allocate_bsorb_pol.f90)
target_link_libraries(bsr_pol bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
