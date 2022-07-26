!======================================================================
      Subroutine receive_results_MPI(process, ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,      only: myid, ierr, noper, joper, JT_oper
      Use coef_list,      only: ntrm, ncoef, idfc, intc, coef
      Use term_exp,       only: kt1, kt2, IP_kt1, IP_kt2
      Use ndet_list,      only: ndet, ldet, KPD, IPD, NPD
      Use ndef_list,      only: ndef, ldef, KPF, IPF, NPF

      Implicit none
      Integer, intent(out) :: process, jc, ic, status(MPI_STATUS_SIZE)

! ... Receive base info
      Call MPI_RECV(ic,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      process = status(MPI_SOURCE)
      Call MPI_RECV(jc,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ncoef,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ntrm,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(joper,noper,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JT_oper)) Deallocate(JT_oper)
      Allocate(JT_oper(ntrm,noper))
      Call MPI_RECV(JT_oper,noper*ntrm,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive term info
      Call MPI_RECV(kt1,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt1)) Deallocate(IP_kt1)
      Allocate(IP_kt1(kt1))
      Call MPI_RECV(IP_kt1,kt1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(kt2,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt2)) Deallocate(IP_kt2)
      Allocate(IP_kt1(kt2))
      Call MPI_RECV(IP_kt2,kt2,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive coefficient info
      if(allocated(idfc)) Deallocate(idfc)
      Allocate(idfc(ncoef))
      Call MPI_RECV(idfc,ncoef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(intc)) Deallocate(intc)
      Allocate(intc(ncoef))
      Call MPI_RECV(intc,ncoef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(coef)) Deallocate(coef)
      Allocate(coef(ntrm,ncoef))
      Call MPI_RECV(coef,ntrm*ncoef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive det info
      Call MPI_RECV(ndet,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldet,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPD)) Deallocate(KPD)
      Allocate(KPD(ndet))
      Call MPI_RECV(KPD,ndet,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPD)) Deallocate(IPD)
      Allocate(IPD(ndet))
      Call MPI_RECV(IPD,ndet,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPD)) Deallocate(NPD)
      Allocate(NPD(ldet))
      Call MPI_RECV(NPD,ldet,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive def info
      Call MPI_RECV(ndef,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldef,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPF)) Deallocate(KPF)
      Allocate(KPF(ndef))
      Call MPI_RECV(KPF,ndef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPF)) Deallocate(IPF)
      Allocate(IPF(ndef))
      Call MPI_RECV(IPF,ndef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPF)) Deallocate(NPF)
      Allocate(NPF(ldef))
      Call MPI_RECV(NPF,ldef,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,ierr)

      End Subroutine receive_results_MPI
