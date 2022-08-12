add_executable(bsr_hd
  mod_hd.f90               bsr_hd.f90              b_out.f90
  diag_hd.f90              h_out.f90               rsol_out.f90
  r_arg.f90                sub1_hd.f90             w_out.f90
  find_channel_label.f90   bsr_hd_inf.f90)
target_link_libraries(bsr_hd bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
