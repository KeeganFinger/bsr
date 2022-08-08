!======================================================================
      Subroutine send_results_MPI(ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,      only: myid, ierr, noper, joper, JT_oper
      Use coef_list,      only: ntrm, ncoef, idfc, intc, coef
      Use term_exp,       only: kt1, kt2, IP_kt1, IP_kt2
      Use ndet_list,      only: ndet, ldet, KPD, IPD, NPD, JPD, idet
      Use ndef_list,      only: ndef, ldef, KPF, IPF, NPF, JPF, idef

      Implicit none
      Integer, intent(in) :: ic, jc

! ... Send base info
      Call MPI_SSEND(ic,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
      if(ic.lt.0) return
      Call MPI_SSEND(jc,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ncoef,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ntrm,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(joper,noper,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(JT_oper,noper*ntrm,MPI_INTEGER,0,5,MPI_COMM_WORLD,ierr)

! ... Send term info
      Call MPI_SSEND(kt1,1,MPI_INTEGER,0,6,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IP_kt1,kt1,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(kt2,1,MPI_INTEGER,0,8,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IP_kt2,kt2,MPI_INTEGER,0,9,MPI_COMM_WORLD,ierr)

! ... Send coefficient info
      Call MPI_SSEND(idfc,ncoef,MPI_INTEGER,0,10,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(intc,ncoef,MPI_INTEGER,0,11,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION,0,12,MPI_COMM_WORLD,ierr)

! ... Send det info
      Call MPI_SSEND(ndet,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ldet,1,MPI_INTEGER,0,14,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(KPD,idet,MPI_INTEGER,0,15,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IPD,idet,MPI_INTEGER,0,16,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(JPD,idet,MPI_INTEGER,0,17,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(NPD,idet,MPI_INTEGER,0,18,MPI_COMM_WORLD,ierr)

! ... Send def info
      Call MPI_SSEND(ndef,1,MPI_INTEGER,0,19,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ldef,1,MPI_INTEGER,0,20,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(KPF,idef,MPI_INTEGER,0,21,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IPF,idef,MPI_INTEGER,0,22,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(JPF,idef,MPI_INTEGER,0,23,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(NPF,idef,MPI_INTEGER,0,24,MPI_COMM_WORLD,ierr)

!      ndet = 0
!      ldet = 0
!      ndef = 0
!      ldef = 0
!      IPF = 0;

      End Subroutine send_results_MPI
