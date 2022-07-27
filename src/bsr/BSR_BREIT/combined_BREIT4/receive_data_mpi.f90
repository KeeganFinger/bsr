!======================================================================
      Subroutine receive_data_MPI(ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,     only: ierr, noper, joper, JT_oper, CT_oper
      Use conf_LS,       only: ne
      Use coef_list,     only: ntrm
      Use spin_orbitals, only: NNsym1, NNsym2, Lsym1, Lsym2
      Use term_exp,      only: kt1, kdt1, ILT1, IST1, MLT, MST, &
                               kt2, kdt2, ILT2, IST2, &
                               IM_det1, IS_det1, &
                               IM_det2, IS_det2, &
                               IP_kt1, IP_kt2, &
                               C_det1, C_det2

      Implicit none
      Integer, intent(out) :: ic, jc
      Integer :: status(MPI_STATUS_SIZE)

! ... Receive base info
      Call MPI_RECV(ic,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)
      if(ic.le.0) Return ! ic controls looping/kill processes

      Call MPI_RECV(jc,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ntrm,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(joper,noper,MPI_INTEGER,0,3,MPI_COMM_WORLD, status,ierr)
      if(allocated(JT_oper)) Deallocate(JT_oper)
      Allocate(JT_oper(ntrm,noper))
      Call MPI_RECV(JT_oper,ntrm*noper,MPI_INTEGER,0,4,MPI_COMM_WORLD, status,ierr)
      if(allocated(CT_oper)) Deallocate(CT_oper)
      Allocate(CT_oper(ntrm, noper))

! ... Receive outer loop info (suffix 1)
      Call MPI_RECV(kt1,1,MPI_INTEGER,0,5,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(kdt1,1,MPI_INTEGER,0,6,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ILT1,1,MPI_INTEGER,0,7,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(IST1,1,MPI_INTEGER,0,8,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(MLT,1,MPI_INTEGER,0,9,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(MST,1,MPI_INTEGER,0,10,MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt1)) Deallocate(IP_kt1)
      Allocate(IP_kt1(kt1))
      Call MPI_RECV(IP_kt1,kt1,MPI_INTEGER,0,11,MPI_COMM_WORLD,status,ierr)
      if(allocated(IM_det1)) Deallocate(IM_det1)
      Allocate(IM_det1(ne,kdt1))
      Call MPI_RECV(IM_det1,ne*kdt1,MPI_INTEGER,0,12,MPI_COMM_WORLD,status,ierr)
      if(allocated(IS_det1)) Deallocate(IS_det1)
      Allocate(IS_det1(ne,kdt1))
      Call MPI_RECV(IS_det1,ne*kdt1,MPI_INTEGER,0,13,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(nnsym1,ne,MPI_INTEGER,0,14,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(lsym1,ne,MPI_INTEGER,0,15,MPI_COMM_WORLD,status,ierr)
      if(allocated(C_det1)) Deallocate(C_det1)
      Allocate(C_det1(kt1,kdt1))
      Call MPI_RECV(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,0,16,MPI_COMM_WORLD,status,ierr)

! ... Receive inner loop info (suffix 2)
      Call MPI_RECV(kt2,1,MPI_INTEGER,0,17,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(kdt2,1,MPI_INTEGER,0,18,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(ILT2,1,MPI_INTEGER,0,19,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(IST2,1,MPI_INTEGER,0,20,MPI_COMM_WORLD,status,ierr)
      if(allocated(IP_kt2)) Deallocate(IP_kt2)
      Allocate(IP_kt2(kt2))
      Call MPI_RECV(IP_kt2,kt2,MPI_INTEGER,0,21,MPI_COMM_WORLD,status,ierr)
      if(allocated(IM_det2)) Deallocate(IM_det2)
      Allocate(IM_det2(ne,kdt2))
      Call MPI_RECV(IM_det2,ne*kdt2,MPI_INTEGER,0,22,MPI_COMM_WORLD,status,ierr)
      if(allocated(IS_det2)) Deallocate(IS_det2)
      Allocate(IS_det2(ne,kdt2))
      Call MPI_RECV(IS_det2,ne*kdt2,MPI_INTEGER,0,23,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(nnsym2,ne,MPI_INTEGER,0,24,MPI_COMM_WORLD,status,ierr)
      Call MPI_RECV(lsym2,ne,MPI_INTEGER,0,25,MPI_COMM_WORLD,status,ierr)
      if(allocated(C_det2)) Deallocate(C_det2)
      Allocate(C_det2(kt2,kdt2))
      Call MPI_RECV(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,0,26,MPI_COMM_WORLD,status,ierr)

      End Subroutine receive_data_MPI
