!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations
!-----------------------------------------------------------------------
      Use bsr_breit
      Use spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      Use term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, MLT,MST, ic_case
      Use conf_LS,       only: ne
      Use symc_list_LS,  only: JC_need, IC_need, nsymc
      Use coef_list,     only: ntrm,ctrm, ncoef, coef, idfc, intc,&
                               ijhm, ctrm, ipcoef
      Use zoef_list,     only: nzoef
      Use ndet_list,     only: KPD,IPD,NPD,JPD,ndet,ldet
      Use ndef_list,     only: KPF,IPF,NPF,JPF,ndef,ldef

      Implicit none
      Integer :: k1,k2, it,jt, MLT2,MST2, i,m,k, is,js,ic,jc
      Real(8) :: C_ee,C_so,C_ss, zero=0.d0, one=1.d0
      Integer, external :: IDEF_cme
      Real(8), external :: Z_3j
      Integer(8) :: ij
      Integer(8), external :: DEF_ij8
      Character(80) :: conf, filename

!----------------------------------------------------------------------
! ... cycle 1 over configurations:

      rewind(nud)
      Do is=1,ic_case

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST

       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

       if(Allocated(IM_det1)) Deallocate(IM_det1)
       Allocate(IM_det1(ne,kdt1)); Read(nud) IM_det1

       if(Allocated(IS_det1)) Deallocate(IS_det1)
       Allocate(IS_det1(ne,kdt1)); Read(nud) IS_det1

       read(nud) NNsym1(1:ne)
       read(nud) Lsym1(1:ne)

       if(IC_need(ic).eq.0) Cycle

       Call CPU_TIME(t1)

       Call Alloc_boef(-1)
       Call Alloc_blk (-1)

!----------------------------------------------------------------------
! ... cycle 2 over configurations:

      rewind(nud)
      Do js=1,is

       Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       if(Allocated(IM_det2)) Deallocate(IM_det2)
       Allocate(IM_det2(ne,kdt2)); Read(nud) IM_det2

       if(Allocated(IS_det2)) Deallocate(IS_det2)
       Allocate(IS_det2(ne,kdt2)); Read(nud) IS_det2

       read(nud) NNsym2(1:ne)
       read(nud) Lsym2(1:ne)

       if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle
       if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle
       ij=DEF_ij8(ic,jc);  if(JC_need(ij).eq.0) Cycle


!----------------------------------------------------------------------
! ...  define number of terms:

       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1)
       Do k2=1,kt2; jt=IP_kt2(k2)
        if(is.eq.js.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do

!----------------------------------------------------------------------
! ...  joper and JT_oper:

       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))
       if(IDEF_cme(is,js).eq.0) Cycle

      ! write(filename,'(a,i3.3,i3.3)') 'data.',is,js
      ! open(100,file=filename)
      !
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

!----------------------------------------------------------------------
! ... define the normalization constants for different operators:

       C_so = zero
       if(joper(4)+joper(5).ne.0) &
        C_so = Z_3j(ILT1,-MLT+2,3,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,3,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
       if(C_so.ne.zero) C_so=one/C_so

       C_ss = zero
       if(joper(6).ne.0) &
        C_ss = Z_3j(ILT1,-MLT+2,5,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,5,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
       if(C_ss.ne.zero) C_ss=one/C_ss

       C_ee = one; if(ILT1.ne.ILT2.or.IST1.ne.IST2) C_ee = zero

       if( abs(C_ee) + abs(C_so) + abs(C_ss) .eq. zero) Cycle

       coper = zero
       if(joper(1).gt.0) coper(1) = C_ee
       if(joper(2).gt.0) coper(2) = C_ee
       if(joper(3).gt.0) coper(3) = C_ee
       if(joper(4).gt.0) coper(4) = C_so
       if(joper(5).gt.0) coper(5) = C_so
       if(joper(6).gt.0) coper(6) = C_ss
       if(joper(7).gt.0) coper(7) = C_ee

!----------------------------------------------------------------------
! ...  initial allocations:

       Call Alloc_coef(-1)

!----------------------------------------------------------------------
! ...  calculations:

       Do kd1 = 1,kdt1

        Msym1(1:ne)=IM_det1(1:ne,kd1)
        Ssym1(1:ne)=IS_det1(1:ne,kd1)

        Call Det_orbitals1

       Do kd2 = 1,kdt2

!        if(is.eq.js.and.kd2.lt.kd1) Cycle       ???

        Msym2(1:ne)=IM_det2(1:ne,kd2)
        Ssym2(1:ne)=IS_det2(1:ne,kd2)

        nzoef = 0;      Call Det_orbitals2
        if(nzoef.gt.0) then
          Call Term_loop(is,js)
        endif

       End do

       if(mod(kd1,1000).eq.0) then
        Call CPU_TIME(t2)
        if(time_limit.gt.0.d0.and.(t2-t0)/60.gt.time_limit) go to 10
       end if

       End do

! ...  store results for given config.s:

      write(filename,'(a,I3.3,I3.3)') 'IPD.',is,js
      open(1001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPD.',is,js
      open(1002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPD.',is,js
      open(1003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPD.',is,js
      open(1004,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'IPF_send.',is,js
      open(2001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPF.',is,js
      open(2002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPF.',is,js
      open(2003,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPF.',is,js
      open(2004,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'idfc.',is,js
!      open(3001,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'intc.',is,js
!      open(3002,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'coef.',is,js
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

!      print *, 'ic,jc,ndet,ldet,ndef,ldef',is,js,ndet,ldet,ndef,ldef
      write(1001,*) IPD(1:ndet)
      write(1002,*) JPD(1:ndet)
      write(1003,*) KPD(1:ndet)
      write(1004,*) NPD(1:ldet)
      write(2001,*) shape(IPF)
      write(2001,*) IPF
      write(2002,*) JPF(1:ndef)
      write(2003,*) KPF(1:ndef)
      write(2004,*) NPF(1:ldef)
!      write(3001,*) idfc(1:ncoef)
!      write(3002,*) intc(1:ncoef)
!      write(3003,*) coef(1:ntrm,1:ncoef)
      write(4001,*) ijhm(1:ntrm)
      write(4002,*) ctrm(1:ntrm)
      write(4003,*) ipcoef(1:ncoef)
      write(5001,*) IP_kt1(1:kt1)
      write(5002,*) IP_kt2(1:kt2)

       Call Add_res(nur,is,js); Call Add_it_oper(is,js)

      End do    ! over js

      Call CPU_TIME(t2)

      Call Symc_conf(ic,conf)

      write(pri,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
        ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec.',trim(conf)
      write(*  ,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
        ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec.',trim(conf)

       if(time_limit.gt.0.d0.and.(t2-t0)/60.gt.time_limit) Exit

      End do    ! over is

   10 Continue

      End Subroutine Conf_loop
