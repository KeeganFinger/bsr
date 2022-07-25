!=====================================================================
!     PROGRAM   B S R _ B R E I T _ M P I                 version: 4
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient in non-orthogonal mode
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!
!    klsp1,klsp2  - range of partial wave in BSR calculations,
!                   then cfg.001, cfg.002, ..., are input files
!                   (default -> 0, with input file is cfg.inp)
!
!    oper  - character(7), where each position can be 0 or 1,
!            and indicate the operator under consideration:
!            oper(1) - OVERLAPS
!            oper(2) - KINATIC ENERGY
!            oper(3) - TWO-ELECTRON ELECTROSTATIC
!            oper(4) - SPIN ORBIT
!            oper(5) - SPIN-OTHER-ORBIT
!            oper(6) - SPIN-SPIN
!            oper(7) - ORBIT-ORBIT
!            Default -> 1110000 - non-relativistic calculations
!
!    mk    - max.multipole index (default -> 9, see module param_br)
!
!----------------------------------------------------------------------
!
!    example:    1.  bsr_breit
!                2.  bsr_breit klsp1=1 klsp2=5 oper=1111110
!                3.  bsr_breit km=5
!
!----------------------------------------------------------------------
!
!    INPUT FILES:
!
!    cfg.nnn     -  configuration list for partial wave nnn = klsp
!                   (cfg.inp in case klsp = 0, default)
!
!    int_bnk.nnn -  input data bank for angular coefficients
!                   (optional; int_bnk in case klsp = 0)
!
!
!    OUTPUT FILES:
!
!    int_bnk.nnn  - output data bank for angular coefficients
!                   (int_bnk in case klsp = 0)
!
!---------------------------------------------------------------------
      Use MPI

      Use bsr_breit
      USE conf_LS,      only: ne
      Use symc_list_LS, only: nsymc
      Use term_exp,     only: ic_case

      Implicit none
      Integer :: l,mls_max

!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      t0 = MPI_WTIME()
      if(myid.eq.0) then
       open(pri,file=AF_p)
       write(*,'(a,i6)') 'nprocs = ', nprocs
       write(pri,'(a,i6)') 'nprocs = ', nprocs
       Allocate(ip_proc(nprocs)); ip_proc=0
      end if

!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;
      Call br_arg

      do klsp=klsp1,klsp2

        if(myid.eq.0) then
          write(pri,'(/a,i5)') 'Partial wave: ',klsp
          Call open_c_file
          Call open_int_inf
          if(new.eq.1) write(pri,'(/a)') 'It is new calculations '
          if(new.eq.0) write(pri,'(/a)') 'It is continued calculations '
          Call Read_conf
        endif

        Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(icalc.eq.0) cycle

        if(myid.eq.0) Call Read_dets(nub,new)

        if(myid.eq.0) then
          Call def_maxl(l)
          mls_max=4*l+2
        endif

        Call MPI_BCAST(mls_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call Alloc_spin_orbitals(ne,mls_max)

        if(myid.eq.0) then
          Call open_det_exp
          Call open_int_int
        endif

        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(myid.eq.0) Call Conf_loop_mpi
        if(myid.ne.0) Call Conf_calc_mpi

        if(myid.eq.0) Call Record_results

        t2=MPI_WTIME()
        if(myid.eq.0) write(pri,'(/a,F12.2,a)') &
          ' Partial wave:',(t2-t0)/60,' min'

      enddo !klsp

      Call MPI_FINALIZE(l)

      END 
