add_executable(bsr_hdb4
	mod_hd.f90 mod_blacs.f90 bsr_hdb.f90
	add_exp.f90 add_mvc.f90 b_out.f90
	diag_hd.f90 find_channel_label.f90 h_out.f90
	pri_mainpar.f90 r_arg.f90 read_diag.f90
	rsol_out.f90 sub1_hd.f90 transform_mat.f90
	transform_ovl.f90 w_out.f90)
target_link_libraries(bsr_hdb4 bs zconfls zcom ${BLAS_LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} MPI::MPI_Fortran)
