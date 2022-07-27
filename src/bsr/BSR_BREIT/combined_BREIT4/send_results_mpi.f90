!======================================================================
      Subroutine send_results_MPI(ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,      only: myid, ierr, noper, joper, JT_oper
      Use coef_list,      only: ntrm, ncoef, idfc, intc, coef
      Use term_exp,       only: kt1, kt2, IP_kt1, IP_kt2
      Use ndet_list,      only: ndet, ldet, KPD, IPD, NPD
      Use ndef_list,      only: ndef, ldef, KPF, IPF, NPF

      Implicit none
      Integer, intent(in) :: ic, jc

! ... Send base info
      Call MPI_SEND(ic,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
      if(ic.lt.0) return
      Call MPI_SEND(jc,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ncoef,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ntrm,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(joper,noper,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(JT_oper,noper*ntrm,MPI_INTEGER,0,5,MPI_COMM_WORLD,ierr)

! ... Send term info
      Call MPI_SEND(kt1,1,MPI_INTEGER,0,6,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt1,kt1,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kt2,1,MPI_INTEGER,0,8,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt2,kt2,MPI_INTEGER,0,9,MPI_COMM_WORLD,ierr)

! ... Send coefficient info
      Call MPI_SEND(idfc,ncoef,MPI_INTEGER,0,10,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(intc,ncoef,MPI_INTEGER,0,11,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION,0,12,MPI_COMM_WORLD,ierr)

! ... Send det info
      Call MPI_SEND(ndet,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ldet,1,MPI_INTEGER,0,14,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(KPD,ndet,MPI_INTEGER,0,15,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPD,ndet,MPI_INTEGER,0,16,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPD,ldet,MPI_INTEGER,0,17,MPI_COMM_WORLD,ierr)

! ... Send def info
      Call MPI_SEND(ndef,1,MPI_INTEGER,0,18,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ldef,1,MPI_INTEGER,0,19,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(KPF,ndef,MPI_INTEGER,0,20,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPF,ndef,MPI_INTEGER,0,21,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPF,ldef,MPI_INTEGER,0,22,MPI_COMM_WORLD,ierr)

      End Subroutine send_results_MPI
