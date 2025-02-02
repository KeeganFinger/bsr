set(BSR_BMAT_COMMON_SOURCES
  MOD_bsr_breit.f90
  MOD_bsr_mat.f90
  alloc_orb_overlaps.f90
  breit_me.f90
  prepare_det_expn.f90
  read_arg.f90
  rw_is_need.f90
  read_conf.f90
  term_loop.f90
  check_idef.f90
  do_coef.f90
  idef_type.f90
  read_data.f90
  open_br.f90)

add_executable(bsr_bmat
  ${BSR_BMAT_COMMON_SOURCES}
  bsr_bmat.f90
  add_res.f90
  conf_loop.f90)
target_link_libraries(bsr_bmat zconfls shells zcom bs ${BLAS_LAPACK_LIBRARIES})

if(MPI_Fortran_FOUND)
   add_executable(bsr_bmat_mpi
     ${BSR_BMAT_COMMON_SOURCES}
     br_arg.f90
     bsr_bmat_mpi.f90
     add_res_mpi.f90
     conf_loop_mpi.f90
     get_res_mpi.f90
     send_res_mpi.f90
     conf_calc.f90
     get_det_exp.f90
     send_det_exp.f90)
   target_link_libraries(bsr_bmat_mpi zconfls shells zcom bs bsr_mpi ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
endif()
