!=======================================================================
      Subroutine Conf_loop_mpi
!=======================================================================
!     run loop over configurations
!-----------------------------------------------------------------------
      Use bsr_breit
      Use spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      Use term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, MLT,MST, ic_case
      Use conf_LS,       only: ne
      Use symc_list_LS,  only: JC_need, IC_need, nsymc
      Use coef_list,     only: ntrm,ctrm, ncoef
      Use zoef_list,     only: nzoef

      Implicit none
      Integer :: k1,k2, it,jt, MLT2,MST2, i,m,k, is,js,ic,jc
      Real(8) :: C_ee,C_so,C_ss, zero=0.d0, one=1.d0
      Integer, external :: IDEF_cme
      Real(8), external :: Z_3j
      Integer(8) :: ij
      Integer(8), external :: DEF_ij8
      Character(80) :: conf

!----------------------------------------------------------------------
! ... cycle 1 over configurations:
      t1=MPI_WTIME()
      rewind(nud)
      Do is=1,ic_case
        Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST

        if(Allocated(IP_kt1)) Deallocate(IP_kt1)
        Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

        if(Allocated(C_det1)) Deallocate(C_det1)
        Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

        if(Allocated(IM_det1)) Deallocate(IM_det1)
        Allocate(IM_det1(ne,kdt1)); Read(nud) IM_det1

        if(Allocated(IS_det1)) Deallocate(IS_det1)
        Allocate(IS_det1(ne,kdt1)); Read(nud) IS_det1

        read(nud) NNsym1(1:ne)
        read(nud) Lsym1(1:ne)

        if(IC_need(ic).eq.0) Cycle

        Call CPU_TIME(t1)

        Call Alloc_boef(-1)
        Call Alloc_blk(-1)

        rewind(nud)
        Do js=1,is

          Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

          if(Allocated(IP_kt2)) Deallocate(IP_kt2)
          Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

          if(Allocated(C_det2)) Deallocate(C_det2)
          Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

          if(Allocated(IM_det2)) Deallocate(IM_det2)
          Allocate(IM_det2(ne,kdt2)); Read(nud) IM_det2

          if(Allocated(IS_det2)) Deallocate(IS_det2)
          Allocate(IS_det2(ne,kdt2)); Read(nud) IS_det2

          read(nud) NNsym2(1:ne)
          read(nud) Lsym2(1:ne)

          if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle
          if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle
          ij=DEF_ij8(ic,jc);  if(JC_need(ij).eq.0) Cycle

!----------------------------------------------------------------------
! ...  define number of terms:

          ntrm = 0
          Do k1=1,kt1; it=IP_kt1(k1)
            Do k2=1,kt2; jt=IP_kt2(k2)
              if(is.eq.js.and.it.gt.jt) Cycle
              ntrm = ntrm + 1
            End do
          End do

!----------------------------------------------------------------------
! ...  joper and JT_oper:

          if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
          Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))
          if(IDEF_cme(is,js).eq.0) Cycle

!----------------------------------------------------------------------
! ...  calculations:

          Do i=1,nprocs-1
            if(ip_proc(i).ne.0) Cycle
            Call Send_det_exp(i,is,js)
            met = i
            ip_proc(i) = 1
            exit
          End do

          if(met.eq.0) then
            Call Get_res(i,is,js)
            Call Add_res_mpi(nur,is,js)

            Call Add_it_oper_mpi(is,js)

            Call Send_det_exp(i,is,js)
          end if

        End do    ! over jc

        t2=MPI_WTIME()

        write(*,'(a,4i8,2f10.2,a,5x,a)') 'ic,ic_total,kt,kdt', iis,ic_case,kt1,kdt1, &
          (t2-t3)/60, (t2-t0)/60, ' min.', conf_is(iis)

! ... release processes
        Do i=1,nprocs-1
          Call Send_det_exp(i,-1,-1)
        End do

      End do    ! over ic

      End Subroutine Conf_loop_mpi
