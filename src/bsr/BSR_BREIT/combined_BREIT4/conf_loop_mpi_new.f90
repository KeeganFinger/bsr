!=======================================================================
      Subroutine Conf_loop_mpi
!=======================================================================
!     run loop over configurations
!-----------------------------------------------------------------------
      Use MPI
      Use bsr_breit
      Use spin_orbitals, only: NNsym1, Lsym1, &
                               NNsym2, Lsym2
      Use term_exp,      only: kt1, kdt1, ILT1, IST1, MLT, MST, &
                               IP_kt1, C_det1, IM_det1, IS_det1, &
                               kt2, kdt2, ILT2, IST2, &
                               IP_kt1, C_det1, IM_det1, IS_det1
      Use conf_LS,       only: ne
      Use symc_list_LS,  only: IC_need, JC_need
      Use coef_list,     only: ntrm

      Implicit none
      Integer, external :: IDEF_cme
      Real(8), external :: Z_3j
      Integer(8), external :: DEF_ij8
      Integer :: MLT2, MST2
      Integer :: is, js, k1, k2, proc ! iterators
      Integer :: ic, jc ! indexes
      Integer :: send ! MPI communicator


! ... Outer loop over configurations

      rewind(nud)
      do is=1,ic_case

        read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST

        if(Allocated(IP_kt1)) Deallocate(IP_kt1)
        Allocate(IP_kt1(kt1))
        Read(nud) IP_kt1

        if(Allocated(C_det1)) Deallocate(C_det1)
        Allocate(C_det1(kt1,kdt1))
        Read(nud) C_det1

        if(Allocated(IM_det1)) Deallocate(IM_det1)
        Allocate(IM_det1(ne,kdt1))
        Read(nud) IM_det1

        if(Allocated(IS_det1)) Deallocate(IS_det1)
        Allocate(IS_det1(ne,kdt1))
        Read(nud) IS_det1

        read(nud) NNsym1(1:ne)
        read(nud) Lsym1(1:ne)

        if(IC_need(ic).eq.0) Cycle

        t1 = MPI_WTIME()

! ... Inner loop over configurations
        rewind(nud)
        do js=1,is

          Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

          if(Allocated(IP_kt2)) Deallocate(IP_kt2)
          Allocate(IP_kt2(kt2))
          Read(nud) IP_kt2

          if(Allocated(C_det2)) Deallocate(C_det2)
          Allocate(C_det2(kt2,kdt2))
          Read(nud) C_det2

          if(Allocated(IM_det2)) Deallocate(IM_det2)
          Allocate(IM_det2(ne,kdt2))
          Read(nud) IM_det2

          if(Allocated(IS_det2)) Deallocate(IS_det2)
          Allocate(IS_det2(ne,kdt2))
          Read(nud) IS_det2

          read(nud) NNsym2(1:ne)
          read(nud) Lsym2(1:ne)

          if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle
          if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle
          ij=DEF_ij8(ic,jc)
          if(JC_need(ij).eq.0) Cycle

! ... Define number of terms
          ntrm = 0
          Do k1=1,kt1
            it=IP_kt1(k1)
            Do k2=1,kt2
              jt=IP_kt2(k2)
              if(is.eq.js.and.it.gt.jt) Cycle
              ntrm = ntrm + 1
            End do
          End do

! ... Setup operators
          if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
          Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))
          if(IDEF_cme(is,js).eq.0) Cycle

! ... Begin calculations

          send = 0
          do proc=1,nprocs-1
            if(proc_status(proc).ne.0) cycle !skip if process proc is busy
            Call send_data_MPI(proc,ic,jc)
            proc_status(proc) = 1 !set process proc status to busy
            send = 1 !sent data, so don't receive yet
            exit !exit loop if data sent to a process
          enddo

          ! receive data if no processes available
          if(send.eq.0) then
            Call receive_results_MPI(proc,ic,jc)
            Call add_res(nur,is,js)
            Call add_it_oper(is,js)
          endif

        enddo ! end inner configuration loop

        t2 = MPI_WTIME()

        Call Symc_conf(ic, conf)

        write(pri,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
          ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec.',trim(conf)
        write(*  ,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
          ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec.',trim(conf)

      enddo ! end outer configuration loop

      End Subroutine conf_loop_MPI
