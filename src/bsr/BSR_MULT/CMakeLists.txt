set(BSR_MULT_COMMON_SOURCES
  mod_coef.f90		mod_zoef.f90		term_loop.f90
  det_mult.f90 		idet_fact.f90		prepare.f90
)

add_executable(bsr_mult
  ${BSR_MULT_COMMON_SOURCES}
  bsr_mult.f90
  mod_param.f90
  mod_spinorb.f90
  add_res.f90
  conf_loop.f90
  def_ic.f90
  rw.f90
  r_conf.f90
  inf_mult.f90
)
target_link_libraries(bsr_mult zconfls shells zcom)

if(MPI_Fortran_FOUND)
  add_executable(bsr_mult_mpi
    ${BSR_MULT_COMMON_SOURCES}
    bsr_mult_mpi.f90
    mod_param_mpi.f90
    mod_spinorb_mpi.f90
    add_res_mpi.f90
    conf_loop_mpi.f90
    def_ic_mpi.f90
    stop_mpi.f90
    read_arg.f90
    br_arg.f90
    send_det_exp.f90
    send_res_mpi.f90
    get_det_exp.f90
    get_res_mpi
    conf_calc.f90
    r_conf_mpi.f90
  )
  target_link_libraries(bsr_mult_mpi zconfls shells zcom MPI::MPI_Fortran)
endif()
    
