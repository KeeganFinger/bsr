!======================================================================
      Subroutine send_data_MPI(process, ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,     only: ierr, noper, joper, JT_oper
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
      Integer, intent(in) :: process, ic, jc

! ... Send base info
      Call MPI_SEND(ic,1,MPI_INTEGER,process,0,MPI_COMM_WORLD,ierr)
      if(ic.le.0) return ! ic controls looping/kill processes

      Call MPI_SEND(jc,1,MPI_INTEGER,process,1,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ntrm,1,MPI_INTEGER,process,2,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(joper,noper,MPI_INTEGER,process,3,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(JT_oper,ntrm*noper,MPI_INTEGER,process,4,MPI_COMM_WORLD,ierr)

! ... Send outer loop info (suffix 1)
      Call MPI_SEND(kt1,1,MPI_INTEGER,process,5,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt1,1,MPI_INTEGER,process,6,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ILT1,1,MPI_INTEGER,process,7,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IST1,1,MPI_INTEGER,process,8,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(MLT,1,MPI_INTEGER,process,9,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(MST,1,MPI_INTEGER,process,10,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt1,kt1,MPI_INTEGER,process,11,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IM_det1,ne*kdt1,MPI_INTEGER,process,12,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IS_det1,ne*kdt1,MPI_INTEGER,process,13,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(nnsym1,ne,MPI_INTEGER,process,14,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(lsym1,ne,MPI_INTEGER,process,15,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,process,16,MPI_COMM_WORLD,ierr)

! ... Send inner loop info (suffix 2)
      Call MPI_SEND(kt2,1,MPI_INTEGER,process,17,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt2,1,MPI_INTEGER,process,18,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ILT2,1,MPI_INTEGER,process,19,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IST2,1,MPI_INTEGER,process,20,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt2,kt2,MPI_INTEGER,process,21,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IM_det2,ne*kdt2,MPI_INTEGER,process,22,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IS_det2,ne*kdt2,MPI_INTEGER,process,23,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(nnsym2,ne,MPI_INTEGER,process,24,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(lsym2,ne,MPI_INTEGER,process,25,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,process,26,MPI_COMM_WORLD,ierr)


      End Subroutine send_data_MPI
