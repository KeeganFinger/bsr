add_executable(dbsr_breit
  mod_dbsr_breit.f90
  mod_nljm.f90
  mod_term_exp.f90
  mod_coef.f90
  mod_boef_list.f90
  dbsr_breit.f90
  add_res.f90
  conf_loop.f90
  det_exp.f90
  det_me.f90
  read_conf_jj.f90
  term_loop.f90
  inf_dbsr_breit.f90
  me_jj.f90)
target_link_libraries(dbsr_breit zconfjj zcom)

if(MPI_Fortran_FOUND)
  add_executable(dbsr_breit_mpi
    dbsr_breit_mpi.f90
    mod_param.f90
    mod_boef.f90
    mod_coef_mpi.f90
    mod_nljm_mpi.f90
    mod_term_exp_mpi.f90
    add_res_mpi.f90
    br_conf_jj.f90
    conf_calc.f90
    conf_loop_mpi.f90
    det_exp_mpi.f90
    det_me_mpi.f90
    get_res.f90
    open_jj.f90
    read_conf_jj_mpi.f90
    term_loop_mpi.f90
  )
  target_link_libraries(dbsr_breit_mpi zconfjj zcom MPI::MPI_Fortran)
endif()
