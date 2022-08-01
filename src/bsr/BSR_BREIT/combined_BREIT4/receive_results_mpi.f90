!======================================================================
      Subroutine receive_results_MPI(process, ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,      only: myid, ierr, noper, joper, JT_oper
      Use coef_list,      only: ntrm, ncoef, idfc, intc, coef, ijhm, &
                                ctrm, ipcoef
      Use term_exp,       only: kt1, kt2, IP_kt1, IP_kt2
      Use ndet_list,      only: ndet, ldet, KPD, IPD, NPD, JPD
      Use ndef_list,      only: ndef, ldef, KPF, IPF, NPF, JPF

      Implicit none
      Integer, intent(out) :: process, ic, jc
      Integer :: status(MPI_STATUS_SIZE)
      Character(80) :: filename

! ... Receive base info
      Call MPI_RECV(ic,1,MPI_INTEGER,MPI_ANY_SOURCE,0,&
                    MPI_COMM_WORLD,status,ierr)
      process = status(MPI_SOURCE)
      Call MPI_RECV(jc,1,MPI_INTEGER,MPI_ANY_SOURCE,1,&
                    MPI_COMM_WORLD,status,ierr)

      write(filename,'(a,I3.3,I3.3)') 'IPD.',ic,jc
      open(1001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPD.',ic,jc
      open(1002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPD.',ic,jc
      open(1003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPD.',ic,jc
      open(1004,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'IPF.',ic,jc
!      open(2001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPF.',ic,jc
      open(2002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPF.',ic,jc
      open(2003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPF.',ic,jc
      open(2004,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'idfc.',ic,jc
!      open(3001,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'intc.',ic,jc
!      open(3002,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'coef.',ic,jc
!      open(3003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'ijhm.',ic,jc       
      open(4001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'ctrm.',ic,jc
      open(4002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'ipcoef.',ic,jc
      open(4003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'IP_kt1.',ic,jc
      open(5001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'IP_kt2.',ic,jc
      open(5002,file=filename)

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



      if(allocated(ijhm)) Deallocate(ijhm)
      Allocate(ijhm(ntrm))
      Call MPI_RECV(ijhm,ntrm,MPI_INTEGER,MPI_ANY_SOURCE,25,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(ctrm)) Deallocate(ctrm)
      Allocate(ctrm(ntrm))
      Call MPI_RECV(ctrm,ntrm,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,26,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(ipcoef)) Deallocate(ipcoef)
      Allocate(ipcoef(ncoef))
      Call MPI_RECV(ipcoef,ncoef,MPI_INTEGER,MPI_ANY_SOURCE,27,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive det info
      Call MPI_RECV(ndet,1,MPI_INTEGER,MPI_ANY_SOURCE,13,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldet,1,MPI_INTEGER,MPI_ANY_SOURCE,14,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPD)) Deallocate(KPD)
      Allocate(KPD(ndet))
      Call MPI_RECV(KPD,ndet,MPI_INTEGER,MPI_ANY_SOURCE,15,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPD)) Deallocate(IPD)
      Allocate(IPD(ndet))
      Call MPI_RECV(IPD,ndet,MPI_INTEGER,MPI_ANY_SOURCE,16,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JPD)) Deallocate(JPD)
      Allocate(JPD(ndet))
      Call MPI_RECV(JPD,ndet,MPI_INTEGER,MPI_ANY_SOURCE,17,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPD)) Deallocate(NPD)
      Allocate(NPD(ldet))
      Call MPI_RECV(NPD,ldet,MPI_INTEGER,MPI_ANY_SOURCE,18,&
                    MPI_COMM_WORLD,status,ierr)

! ... Receive def info
      Call MPI_RECV(ndef,1,MPI_INTEGER,MPI_ANY_SOURCE,19,&
                    MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ldef,1,MPI_INTEGER,MPI_ANY_SOURCE,20,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(KPF)) Deallocate(KPF)
      Allocate(KPF(ndef))
      Call MPI_RECV(KPF,ndef,MPI_INTEGER,MPI_ANY_SOURCE,21,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(IPF)) Deallocate(IPF)
      Allocate(IPF(ndef))
      Call MPI_RECV(IPF,ndef,MPI_INTEGER,MPI_ANY_SOURCE,22,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(JPF)) Deallocate(JPF)
      Allocate(JPF(ndef))
      Call MPI_RECV(JPF,ndef,MPI_INTEGER,MPI_ANY_SOURCE,23,&
                    MPI_COMM_WORLD,status,ierr)
      if(allocated(NPF)) Deallocate(NPF)
      Allocate(NPF(ldef))
      Call MPI_RECV(NPF,ldef,MPI_INTEGER,MPI_ANY_SOURCE,24,&
                    MPI_COMM_WORLD,status,ierr)

!      print *, 'ic,jc,ndet,ldet,ndef,ldef',ic,jc,ndet,ldet,ndef,ldef
      write(1001,*) IPD
      write(1002,*) JPD
      write(1003,*) KPD
      write(1004,*) NPD
!      write(2001,*) IPF
      write(2002,*) JPF
      write(2003,*) KPF
      write(2004,*) NPF
!      write(3001,*) idfc
!      write(3002,*) intc
!      write(3003,*) coef
      write(4001,*) ijhm
      write(4002,*) ctrm
      write(4003,*) ipcoef
      write(5001,*) IP_kt1
      write(5002,*) IP_kt2

      End Subroutine receive_results_MPI
