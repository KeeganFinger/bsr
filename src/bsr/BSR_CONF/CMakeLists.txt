set(BSR_CONF_COMMON_SOURCES
  mod_bsr_conf.f90
  check_cfg.f90
  check_pertuber.f90
  def_orth_cond.f90
  def_phys_targ.f90
  gen_conf.f90
  make_coupling.f90
  record_orth.f90
  sub_jj.f90
  sub_jk.f90
  sub_ls.f90
  read_arg.f90
  bsr_conf_inf.f90)

#MESSAGE(WARNING "bsr_conf3: check_cfg.f90 does not compile with libs_04")

add_executable(bsr_conf
  ${BSR_CONF_COMMON_SOURCES}
  bsr_conf.f90)
target_link_libraries(bsr_conf zconfls zcom shells)

#message(WARNING "bsr_conf3_mpi does not build correctly")
if(MPI_Fortran_FOUND)
 add_executable(bsr_conf_mpi
   ${BSR_CONF_COMMON_SOURCES}
   bsr_conf_mpi.f90
   br_arg.f90
   read_orth.f90)
 target_link_libraries(bsr_conf_mpi bsr_mpi zconfls zcom shells MPI::MPI_Fortran)
endif()
