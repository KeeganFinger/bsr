add_executable(bsr_prep
  mod_bsr_prep.f90        bsr_prep.f90            bsr_prep_inf.f90
  check_bsw_file.f90      add_sub_orbital.f90     assign_index.f90
  check_bsr_name.f90      read_arg.f90
  sub_check_orb.f90       sub_phys_orb.f90        sub_phys_pert.f90)
target_link_libraries(bsr_prep bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
