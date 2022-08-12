!=======================================================================
      SUBROUTINE allocate_bsorb_old(m)
!=======================================================================
!     This program allocates (deallocates) space for list of atomic
!     orbitals or reallocates it if necessary
!-----------------------------------------------------------------------
      USE spline_orbitals
      USE spline_param
      USE bsr_pol

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarr1(:), iarr2(:,:)
      Real(8), allocatable :: rarr2(:,:)
      CHARACTER(4), allocatable :: abs(:)

      if(m.le.0) then
      
       if(Allocated(NBS)) &
         Deallocate (NBS,LBS,KBS,MBS,iech,EBS,IBORT,PBS,QBS,OBS)
         nbf = 0;  mbf = 0
      
      elseif(m.gt.mbf.and.nbf.eq.0) then

       if(Allocated(nbs)) & 
          Deallocate (nbs,lbs,kbs,mbs,iech,ebs,IBORT,PBS,QBS,OBS)
       
       mbf = m
       Allocate(nbs(mbf),lbs(mbf),kbs(mbf),ebs(mbf),mbs(1:mbf), &
                iech(1:mbf),IBORT(1:mbf,1:mbf),PBS(1:ns,1:mbf), &
                QBS(1:ns,1:mbf),OBS(1:mbf,1:mbf))
       nbs = 0; lbs = 0; kbs = 0; ebs = '****'; mbs = 0; iech = 0
       IBORT = 2; PBS = 0.d0; QBS = 0.d0; OBS = 0.d0

      elseif(nbf.gt.0.and.m.gt.mbf) then

       Allocate(iarr1(nbf))
       iarr1(1:nbf)=nbs(1:nbf); Deallocate(nbs);  Allocate(nbs(m))
       nbs=0;  nbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=lbs(1:nbf); Deallocate(lbs);  Allocate(lbs(m))
       lbs=0;  lbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=kbs(1:nbf); Deallocate(kbs);  Allocate(kbs(m))
       kbs=0;  kbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=mbs(1:nbf); Deallocate(mbs);  Allocate(mbs(m))
       mbs=0;  mbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=iech(1:nbf); Deallocate(iech);  Allocate(iech(m))
       iech=0; iech(1:nbf)=iarr1(1:nbf)
       Deallocate(iarr1)

       Allocate(abs(nbf))
       abs(1:nbf)=ebs(1:nbf); Deallocate(ebs);  Allocate(ebs(m))
       ebs='****'; ebs(1:nbf)=abs(1:nbf); Deallocate(abs)

       Allocate(iarr2(nbf,nbf))
       iarr2(1:nbf,1:nbf)=ibort(1:nbf,1:nbf); Deallocate(ibort)
       Allocate(ibort(m,m))
       ibort=2; ibort(1:nbf,1:nbf)=iarr2(1:nbf,1:nbf)
       Deallocate(iarr2)

       Allocate(rarr2(ns,nbf))
       rarr2(1:ns,1:nbf)=pbs(1:ns,1:nbf); Deallocate(pbs)
       Allocate(pbs(ns,m))
       pbs=0.d0; pbs(1:ns,1:nbf)=rarr2(1:ns,1:nbf)
       rarr2(1:ns,1:nbf)=qbs(1:ns,1:nbf); Deallocate(qbs)
       Allocate(qbs(ns,m))
       qbs=0.d0; qbs(1:ns,1:nbf)=rarr2(1:ns,1:nbf)
       Deallocate(rarr2)

       Allocate(rarr2(nbf,nbf))
       rarr2(1:nbf,1:nbf)=obs(1:nbf,1:nbf); Deallocate(obs)
       Allocate(obs(m,m))
       obs=0.d0; obs(1:nbf,1:nbf)=rarr2(1:nbf,1:nbf)
       Deallocate(rarr2)

       mbf=m;  ! write(*,*) 'realoc_BS_orb: mbf=',mbf

      end if

      m_borb = 6*mbf + mbf*mbf + 2 * (2*ns*mbf + mbf*mbf) + 6
       
      END SUBROUTINE allocate_bsorb_old
