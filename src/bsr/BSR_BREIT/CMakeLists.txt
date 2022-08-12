set(BSR_BREIT_COMMON_SOURCES
  MOD_bsr_breit.f90         add_res.f90              breit_me.f90
  open_files.f90            read_arg.f90             read_conf.f90
  prepare_it_oper.f90       term_loop.f90
)

add_executable(bsr_breit
  ${BSR_BREIT_COMMON_SOURCES}
  bsr_breit.f90
  conf_loop.f90
)
target_link_libraries(bsr_breit zconfls shells zcom)

if(MPI_Fortran_FOUND)
  add_executable(bsr_breit_mpi
    ${BSR_BREIT_COMMON_SOURCES}
    bsr_breit_mpi_new.f90
    br_arg.f90
    conf_calc_mpi_new.f90
    conf_loop_mpi_new.f90
    send_data_mpi.f90
    send_results_mpi.f90
    receive_data_mpi.f90
    receive_results_mpi.f90
    cache_data.f90
    )
  target_link_libraries(bsr_breit_mpi zconfls shells zcom MPI::MPI_Fortran)
endif()
