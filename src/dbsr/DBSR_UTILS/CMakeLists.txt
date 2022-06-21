#message(WARNING "coef_jj does not link correctly")
add_executable(coef_jj coef_jj.f90)
target_link_libraries(coef_jj zconfjj zcom)

add_executable(dbound_bsw dbound_bsw.f90)
target_link_libraries(dbound_bsw dbs zconfjj zcom)

add_executable(dbound_tab dbound_tab.f90)
target_link_libraries(dbound_tab zcom)

add_executable(dbsr_merge dbsr_merge.f90)
target_link_libraries(dbsr_merge dbs zconfjj zcom)

add_executable(dbsw_tab dbsw_tab.f90)
target_link_libraries(dbsw_tab dbs zconfjj zcom)

add_executable(mult_jj mult_jj.f90)
target_link_libraries(mult_jj zconfjj zcom)

add_executable(bsw_rw bsw_rw.f90)
target_link_libraries(bsw_rw dbs zconfjj zcom)

add_executable(bsrw123 bsrw123.f90)
target_link_libraries(bsrw123 dbs zconfjj zcom)

add_executable(rw_bsw rw_bsw.f90)
target_link_libraries(rw_bsw dbs zconfjj zcom)

add_executable(rw123 rw123.f90)
target_link_libraries(rw123 zconfjj zcom)

add_executable(rw_dat rw_dat.f90)
target_link_libraries(rw_dat zconfjj zcom)
