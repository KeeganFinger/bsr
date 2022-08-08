!======================================================================
!     PROGRAM       B S R _ M A T              version 3
!
!               C O P Y R I G H T -- 2010
!
!     Written by:   Oleg Zatsarinny
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Generate the interaction matrixes in B-spline representation
!======================================================================
!
!   INPUT FILES:
!
!     target         -  description of target states and channels
!     target.bsw     -  target w.f.'s in B-spline basis
!     target_orb     -  list of target physical orbitals
!     knot.dat       -  B-spline grid
!     bsr_par        -  input parameters
!     cfg.nnn        -  configuration list for partial wave nnn
!     int_int.nnn    -  angular coefficient data bank (+ int_inf.nnn)
!     pert_nnn.bsw   -  perturber orbitals if any
!
!   OUTPUT FILES:
!
!     bsr_mat.log   -  general running information
!     mat_log.nnn   -  running information for partial wave nnn
!     bsr_mat.nnn   -  resulting overlap and interaction matrixes
!
!=====================================================================
      Use MPI
      Use bsr_mat
      Use conf_LS

      Implicit none
      Real(8) :: t1,t2,t3
      Integer :: i

! ... Initialize MPI
      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

! ... General output
      if(myid.eq.0) Call bsr_mat_inf
      if(myid.eq.0) then
        open(pri,file=AF_p)
        write(pri,'(a,i6)') 'nprocs = ', nprocs
      end if

! ... Prepare B-spline
      if(myid.eq.0) Call define_grid(z)
      Call br_grid
      Call define_spline

! ... Target information
      if(myid.eq.0) then
        open(nut,file=AF_tar,status='OLD')
        Call R_target(nut)
      endif
      Call br_target

! ... Read input arguments/parameters
      if(myid.eq.0) then
        open(nup,file=AF_par,status='OLD')
        Call Read_arg(nup)
      endif
      Call br_arg

! ... Find common core energy
      if(myid.eq.0) then
        Open(nuw,file=AF_bsw,STATUS='OLD',form='UNFORMATTED')
        Call R_bwfn(nuw)
        Close(nuw)
        if(nwt.ne.nbf) &
         Call Stop_mpi (0,0,' nwt (in target) <> nbf (target.bsw)')
        write(ALSP,'(i3.3)') klsp1
        i=LEN_TRIM(AF_cfg); ; write(AF_cfg(i-2:i),'(i3.3)') klsp1
        Open(nuc,file=AF_cfg,status='OLD')
        Call R_closed(nuc)
        kclosd=nclosd
        Call Bcore
        write(prj,'(/a,i4,a)') &
          'nclosd  =',nclosd,' - common core shells'
        write(prj,'(/a,F15.8,a)') &
          'Bcore   =', EC,'  -  calculated core energy'
        Call Read_rpar(nup,'Ecore',EC)
        write(prj,'(/a,F15.8,a)') &
          'Ecore   =', EC,'  -  Used core energy'
      endif
      Call br_core

! ... Loop over partial waves
      do klsp = klsp1, klsp2
        if(myid.eq.0) write(*,'(/a,i3/)') 'BSR_MAT:  klsp =', klsp

        t1 = MPI_WTIME()
        Call SUB1_MPI
        t2 = MPI_WTIME()

        if(myid.eq.0) then
          write(pri,'(a,5x,T20,f10.1,a)') &
            'Total time:', (t2-t1)/60, ' min'
          write(*,'(a,5x,T20,f10.1,a)') &
            'Total time:', (t2-t1)/60, ' min'
        endif
      enddo

      Call MPI_FINALIZE(ierr)

      End
