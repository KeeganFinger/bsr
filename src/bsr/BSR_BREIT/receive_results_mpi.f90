!======================================================================
      Subroutine receive_results_MPI(process, ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,      only: myid, ierr, noper, joper, JT_oper
      Use coef_list,      only: ntrm, ncoef, idfc, intc, coef
      Use term_exp,       only: kt1, kt2, IP_kt1, IP_kt2
      Use ndet_list,      only: ndet, ldet, KPD, IPD, NPD, JPD, idet
      Use ndef_list,      only: ndef, ldef, KPF, IPF, NPF, JPF, idef

      Implicit none
      Integer, intent(out) :: process, ic, jc
      Integer :: status(MPI_STATUS_SIZE)

! ... Receive base info
      Call MPI_RECV(ic,1,MPI_INTEGER,MPI_ANY_SOURCE,0,&
                    MPI_COMM_WORLD,status,ierr)
      process = status(MPI_SOURCE)
      Call MPI_RECV(jc,1,MPI_INTEGER,MPI_ANY_SOURCE,1,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ncoef,1,MPI_INTEGER,MPI_ANY_SOURCE,2,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ntrm,1,MPI_INTEGER,MPI_ANY_SOURCE,3,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(joper,noper,MPI_INTEGER,MPI_ANY_SOURCE,4,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JT_oper)) Deallocate(JT_oper)
      Allocate(JT_oper(ntrm,noper))
      Call MPI_RECV(JT_oper,noper*ntrm,MPI_INTEGER,MPI_ANY_SOURCE,5,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive term info
      Call MPI_RECV(kt1,1,MPI_INTEGER,MPI_ANY_SOURCE,6,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt1)) Deallocate(IP_kt1)
      Allocate(IP_kt1(kt1))
      Call MPI_RECV(IP_kt1,kt1,MPI_INTEGER,MPI_ANY_SOURCE,7,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(kt2,1,MPI_INTEGER,MPI_ANY_SOURCE,8,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt2)) Deallocate(IP_kt2)
      Allocate(IP_kt2(kt2))
      Call MPI_RECV(IP_kt2,kt2,MPI_INTEGER,MPI_ANY_SOURCE,9,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive coefficient info
      if(allocated(idfc)) Deallocate(idfc)
      Allocate(idfc(ncoef))
      Call MPI_RECV(idfc,ncoef,MPI_INTEGER,MPI_ANY_SOURCE,10,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(intc)) Deallocate(intc)
      Allocate(intc(ncoef))
      Call MPI_RECV(intc,ncoef,MPI_INTEGER,MPI_ANY_SOURCE,11,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(coef)) Deallocate(coef)
      Allocate(coef(ntrm,ncoef))
      Call MPI_RECV(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,12,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive det info
      Call MPI_RECV(ndet,1,MPI_INTEGER,MPI_ANY_SOURCE,13,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldet,1,MPI_INTEGER,MPI_ANY_SOURCE,14,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPD)) Deallocate(KPD)
      Allocate(KPD(idet))
      Call MPI_RECV(KPD,idet,MPI_INTEGER,MPI_ANY_SOURCE,15,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPD)) Deallocate(IPD)
      Allocate(IPD(idet))
      Call MPI_RECV(IPD,idet,MPI_INTEGER,MPI_ANY_SOURCE,16,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JPD)) Deallocate(JPD)
      Allocate(JPD(idet))
      Call MPI_RECV(JPD,idet,MPI_INTEGER,MPI_ANY_SOURCE,17,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPD)) Deallocate(NPD)
      Allocate(NPD(idet))
      Call MPI_RECV(NPD,idet,MPI_INTEGER,MPI_ANY_SOURCE,18,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive def info
      Call MPI_RECV(ndef,1,MPI_INTEGER,MPI_ANY_SOURCE,19,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldef,1,MPI_INTEGER,MPI_ANY_SOURCE,20,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPF)) Deallocate(KPF)
      Allocate(KPF(idef))
      Call MPI_RECV(KPF,idef,MPI_INTEGER,MPI_ANY_SOURCE,21,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPF)) Deallocate(IPF)
      Allocate(IPF(idef))
      Call MPI_RECV(IPF,idef,MPI_INTEGER,MPI_ANY_SOURCE,22,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JPF)) Deallocate(JPF)
      Allocate(JPF(idef))
      Call MPI_RECV(JPF,idef,MPI_INTEGER,MPI_ANY_SOURCE,23,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPF)) Deallocate(NPF)
      Allocate(NPF(idef))
      Call MPI_RECV(NPF,idef,MPI_INTEGER,MPI_ANY_SOURCE,24,&
                    MPI_COMM_WORLD,status,ierr)

      End Subroutine receive_results_MPI
