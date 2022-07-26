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
      Use conf_LS,      only: ne
      Use det_list,     only: ndet,ldet
      Use def_list,     only: ndef,ldef

      Implicit none
      Integer :: i,ii,l,ml
      Real(8) :: tt1,tt2, ttt, total_time

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(myid.eq.0) then
        open(pri,file=AF_p)
        write(*,'(a,i6)') 'nprocs = ', nprocs
        write(pri,'(a,i6)') 'nprocs = ', nprocs
        Allocate(proc_status(nprocs)); proc_status=0
      end if

! ... Read and broadcast command line arguments
      if(myid.eq.0) then
        Open(pri,file=AF_p)
        write(pri,'(/20x,a/20x,a/20x,a/)') &
               '=====================================',     &
               ' CALCULATION OF ANGULAR COEFFICIENTS ',     &
               '====================================='
        Call Read_arg
      endif
      Call br_arg

      total_time = 0.0d0
      Do klsp = klsp1, klsp2
        tt1 = MPI_WTIME()

        if(myid.eq.0) then
          write(pri,'(80(''-''))')
          write(pri,'(/a,i5)') 'Partial wave: ',klsp
          (*  ,'(/a,i5)') 'Partial wave: ',klsp
          Call open_c_file
          Call open_int_inf
          if(new.eq.1) then
            write(pri,'(/a)') 'It is new calculations'
          else
            'It is continued calculations'
          endif
        endif

! ... Read configuration list
        if(myid.eq.0) Call read_conf

        Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(icalc.eq.0) cycle

! ... Extract old results
        if(myid.eq.0) then
          Call Read_dets(nub,new)
          Call def_maxl(l)
          mls_max=4*l+s
        endif

        Call MPI_BCAST(mls_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        Call Alloc_spin_orbitals(ne,mls_max)

! ... Prepare det expansions

        if(myid.eq.0) then
          Call open_det_exp
          Call open_int_int
        endif

! ... Calculate angular symmetry coefficients
        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) Call conf_loop_MPI
        if(myid.ne.0) Call conf_calc_MPI

! ... Record results and runtime data
        if(myid.eq.0) Call record_results

        t2=MPI_WTIME()
        if(myid.eq.0) write(pri,'(/a,F12.2,a)') &
          'Partial wave:',(t2-t1)/60,' min'

      enddo
      Call MPI_FINALIZE()
      End

!======================================================================
      Subroutine Read_dets(nub,new)
!======================================================================
      Implicit none

      Integer :: nub, new

      if(new.eq.1) then
       Call Alloc_det(-1)
       Call Alloc_def(-1)
      else
       Call Load_det(nub)
       Call Load_def(nub)
      end if

      End Subroutine Read_dets


!======================================================================
      Subroutine Record_results
!======================================================================
      Use bsr_breit
      Use det_list
      Use def_list

      Implicit none

      rewind(nub)
      Call Write_symc_LS(nub)
      Call Write_symt_LS(nub)
      Call Record_oper_LS(nub)

      adet=ldet;
      if(ndet.gt.0) adet=adet/ndet
      Call Write_det(nub)
      adef=ldef;
      if(ndef.gt.0) adef=adef/ndef
      Call Write_def(nub)

      close(nub)

! ... print the main dimensions:

      write(pri,'(/a/)') &
          ' Results for new angular symmetry calculations:'
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap determinants =', ndet,adet,ldet
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap factors      =', ndef,adef,ldef
      write(pri,'(a,i10)') &
          ' new coeff.s                    =', nc_new


      End Subroutine Record_results


!======================================================================
      Subroutine Debug_printing
!======================================================================
      Use MPI

      Use bsr_breit
      Use coef_list
      Use zoef_list
      Use boef_list

      Integer :: nn, mt, nreloc, mreloc
      Real(8) :: mem, meb

      if(pri.gt.0) write(pri,'(/a/)') 'debug printing:'

      Call MPI_REDUCE(mem_max_zoef,mem,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_zoef,nn,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(zoef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_zoef = ', mem
       write(pri,'(a,2i10)')  'max_zoef = ', nn, izoef
       write(pri,'(a,i10/)')  'max_zoef = ', nreloc
      end if

      Call MPI_REDUCE(mem_max_coef,mem,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_coef,    nn,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_term,    nt,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(coef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_coef = ', mem
       write(pri,'(a,2i10)')  'max_coef = ', nn, icoef
       write(pri,'(a,2i10)')  'max_term = ', nt
       write(pri,'(a,i10/)')  'realloc  = ', nreloc
      end if

      Call MPI_REDUCE(mem_max_boef,mem,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(mem_max_blk, meb,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_boef,    nn,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_blk,     nt,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(boef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(blk_realloc, mreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_boef = ', mem
       write(pri,'(a,f10.1)') 'mem_blk  = ', meb
       write(pri,'(a,2i10)')  'max_boef = ', nn, iboef
       write(pri,'(a,2i10)')  'max_blk  = ', nt
       write(pri,'(a,i10/)')  'al_boef  = ', nreloc
       write(pri,'(a,i10/)')  'al_blk   = ', mreloc
      end if

      End Subroutine Debug_printing
