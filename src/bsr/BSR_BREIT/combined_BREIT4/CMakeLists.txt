set(BSR_BREIT_COMMON_SOURCES
  MOD_bsr_breit.f90         breit_me.f90
  open_files.f90            read_arg.f90             read_conf.f90
  prepare_it_oper.f90       term_loop.f90
)

add_executable(bsr_breit4
  ${BSR_BREIT_COMMON_SOURCES}
  bsr_breit.f90
  conf_loop.f90
  add_res.f90
)
target_link_libraries(bsr_breit4 zconfls shells zcom)

if(MPI_Fortran_FOUND)
  add_executable(bsr_breit4_mpi
    ${BSR_BREIT_COMMON_SOURCES}
    add_res_mpi.f90
    add_it_oper_mpi.f90
    bsr_breit_mpi_new.f90
    br_arg.f90
    conf_calc_mpi_new.f90
    conf_loop_mpi_new.f90
    get_det_mpi.f90
    get_res_mpi.f90
    send_det_mpi.f90
    send_res_mpi.f90
    )
  target_link_libraries(bsr_breit4_mpi zconfls shells zcom MPI::MPI_Fortran)
endif()
