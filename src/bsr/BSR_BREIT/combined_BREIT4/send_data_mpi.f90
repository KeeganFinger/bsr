!======================================================================
      Subroutine send_data_MPI(process, ic, jc)
!======================================================================
      Use MPI
      Use bsr_breit,     only: ierr, noper, joper, JT_oper
      Use conf_LS,       only: ne
      Use coef_list,     only: ntrm
      Use spin_orbitals, only: NNsym1, NNsym2, Lsym1, Lsym2
      Use ndet_list,     only: ndet, ldet, IPD, JPD, KPD, NPD
      Use ndef_list,     only: ndef, ldef, IPF, JPF, KPF, NPF
      Use term_exp,      only: kt1, kdt1, ILT1, IST1, MLT, MST, &
                               kt2, kdt2, ILT2, IST2, &
                               IM_det1, IS_det1, &
                               IM_det2, IS_det2, &
                               IP_kt1, IP_kt2, &
                               C_det1, C_det2

      Implicit none
      Integer, intent(in) :: process, ic, jc
      ! Character(80) :: filename
      ! write(filename,'(a,i3.3,i3.3)') 'data.',ic,jc
      ! open(100,file=filename)

! ... Send base info
      Call MPI_SSEND(ic,1,MPI_INTEGER,process,0,MPI_COMM_WORLD,ierr)
      if(ic.le.0) return ! ic controls looping/kill processes

      Call MPI_SSEND(jc,1,MPI_INTEGER,process,1,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ntrm,1,MPI_INTEGER,process,2,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(joper,noper,MPI_INTEGER,process,3,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(JT_oper,ntrm*noper,MPI_INTEGER,process,4,MPI_COMM_WORLD,ierr)

! ... Send outer loop info (suffix 1)
      Call MPI_SSEND(kt1,1,MPI_INTEGER,process,5,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(kdt1,1,MPI_INTEGER,process,6,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ILT1,1,MPI_INTEGER,process,7,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IST1,1,MPI_INTEGER,process,8,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(MLT,1,MPI_INTEGER,process,9,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(MST,1,MPI_INTEGER,process,10,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IP_kt1,kt1,MPI_INTEGER,process,11,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IM_det1,ne*kdt1,MPI_INTEGER,process,12,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IS_det1,ne*kdt1,MPI_INTEGER,process,13,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(nnsym1,ne,MPI_INTEGER,process,14,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(lsym1,ne,MPI_INTEGER,process,15,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,process,16,MPI_COMM_WORLD,ierr)

! ... Send inner loop info (suffix 2)
      Call MPI_SSEND(kt2,1,MPI_INTEGER,process,17,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(kdt2,1,MPI_INTEGER,process,18,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(ILT2,1,MPI_INTEGER,process,19,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IST2,1,MPI_INTEGER,process,20,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IP_kt2,kt2,MPI_INTEGER,process,21,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IM_det2,ne*kdt2,MPI_INTEGER,process,22,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(IS_det2,ne*kdt2,MPI_INTEGER,process,23,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(nnsym2,ne,MPI_INTEGER,process,24,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(lsym2,ne,MPI_INTEGER,process,25,MPI_COMM_WORLD,ierr)
      Call MPI_SSEND(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,process,26,MPI_COMM_WORLD,ierr)

      ! write(100,*) 'Base Info:'
      ! write(100,*) ic,jc,ntrm
      ! write(100,*) joper
      ! write(100,*) JT_oper
      ! write(100,*) 'Outer Loop Info:'
      ! write(100,*) kt1,kdt1,ILT1,IST1
      ! write(100,*) MLT,MST
      ! write(100,*) IP_kt1
      ! write(100,*) IM_det1
      ! write(100,*) IS_det1
      ! write(100,*) nnsym1
      ! write(100,*) lsym1
      ! write(100,*) C_det1
      ! write(100,*) 'Inner Loop Info:'
      ! write(100,*) kt2,kdt2,ILT2,IST2
      ! write(100,*) IP_kt2
      ! write(100,*) IM_det2
      ! write(100,*) IS_det2
      ! write(100,*) nnsym2
      ! write(100,*) lsym2
      ! write(100,*) C_det2

! ... Send determinant info
!      Call MPI_SSEND(ndet,1,MPI_INTEGER,process,27,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(ldet,1,MPI_INTEGER,process,28,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(KPD,ndet,MPI_INTEGER,process,29,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(IPD,ndet,MPI_INTEGER,process,30,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(JPD,ndet,MPI_INTEGER,process,31,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(NPD,ldet,MPI_INTEGER,process,32,MPI_COMM_WORLD,ierr)

!      Call MPI_SSEND(ndef,1,MPI_INTEGER,process,33,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(ldef,1,MPI_INTEGER,process,34,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(KPF,ndef,MPI_INTEGER,process,35,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(IPF,ndef,MPI_INTEGER,process,36,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(JPF,ndef,MPI_INTEGER,process,37,MPI_COMM_WORLD,ierr)
!      Call MPI_SSEND(NPF,ldef,MPI_INTEGER,process,38,MPI_COMM_WORLD,ierr)

      End Subroutine send_data_MPI
