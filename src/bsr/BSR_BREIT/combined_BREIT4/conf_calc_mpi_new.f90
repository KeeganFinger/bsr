!=======================================================================
      Subroutine Conf_calc_mpi
!=======================================================================
      Use MPI
      Use bsr_breit
      Use spin_orbitals, only: Msym1, Ssym1, &
                               Msym2, Ssym2
      Use term_exp,      only: kdt1, ILT1, IST1, MLT, MST, &
                               kdt2, ILT2, IST2, &
                               IM_det1, IS_det1, &
                               IM_det2, IS_det2, kt1
      Use conf_LS,       only: ne
      Use zoef_list,     only: nzoef

      Implicit none
      Real(8), parameter :: zero=0.d0, one=1.d0
      Real(8) :: C_ee, C_so, C_ss !normalization factors
      Real(8), external :: Z_3j
      Integer :: kd1, kd2 !iterators
      Integer :: is, js !indexes

! ... Prepare to receive data
      Call Alloc_boef(-1)
      Call Alloc_blk (-1)

1     Call receive_data_MPI(is,js)
      if(is.lt.0) return
      print *, myid,'conf_calc', is, js, kt1

! ... Define normalization constants for different operators
      C_so = zero
      if(joper(4)+joper(5).eq.0) then
        C_so = Z_3j(ILT1,-MLT+2,3,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,3,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
      endif
      if(C_so.ne.zero) C_so = one/C_so

      C_ss = zero
      if(joper(6).ne.0) then
        C_ss = Z_3j(ILT1,-MLT+2,5,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,5,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
      endif
      if(C_ss.ne.zero) C_ss = one/C_ss

      C_ee = one
      if(ILT1.ne.ILT2.or.IST1.ne.IST2) C_ee = zero

      if(abs(C_ee)+abs(C_so)+abs(C_ss).eq.zero) then
        Call send_results_mpi(is,js)
        go to 1
      endif


      coper = zero
      if(joper(1).gt.0) coper(1) = C_ee
      if(joper(2).gt.0) coper(2) = C_ee
      if(joper(3).gt.0) coper(3) = C_ee
      if(joper(4).gt.0) coper(4) = C_so
      if(joper(5).gt.0) coper(5) = C_so
      if(joper(6).gt.0) coper(6) = C_ss
      if(joper(7).gt.0) coper(7) = C_ee

! ... Allocation initialization
      Call Alloc_coef(-1)

! ... Begin calculations
      do kd1 = 1,kdt1
        Msym1(1:ne) = IM_det1(1:ne,kd1)
        Ssym1(1:ne) = IS_det1(1:ne,kd1)

        Call Det_orbitals1

        Do kd2 = 1,kdt2
          Msym2(1:ne) = IM_det2(1:ne,kd2)
          Ssym2(1:ne) = IS_det2(1:ne,kd2)

          nzoef = 0
          Call Det_orbitals2
          if(nzoef.gt.0) Call Term_loop(is,js)
          
        enddo ! loop over kd2

        t3 = MPI_WTIME()
        if(time_limit.gt.0.d0.and.(t3-t0)/60.gt.time_limit*1.1) then
          Call send_results_MPI(-1*is,js)
          go to 1
        endif
      enddo ! loop over kd1

      print *, myid,'conf_calc2', is, js, kt1
      Call send_results_MPI(is,js)
      go to 1

      End Subroutine Conf_calc_MPI
