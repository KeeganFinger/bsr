add_executable(dbsr_mult
  mod_dbsr_mult.f90          mod_nljm.f90               mod_term_exp.f90
  mod_coef.f90               check_mult_bnk.f90         add_res.f90
  conf_loop.f90              dbsr_mult.f90              det_exp.f90
  det_mult.f90               read_arg.f90               term_loop.f90
  Idet_fact.f90              me_jj.f90)
target_link_libraries(dbsr_mult zconfjj zcom)
