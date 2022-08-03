!======================================================================
      Subroutine Add_res(nu,is,js)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      Use bsr_breit
      Use coef_list
      Use term_exp,  only: kt1,kt2, IP_kt1,IP_kt2
      Use ndef_list, only: IPF,JPF,KPF,NPF,ndef,ldef
      Use ndet_list, only: IPD,JPD,KPD,NPD,ndet,ldet
      Use spin_orbitals, only: Lsym1, Msym1, Ssym1,&
                               Lsym2, Msym2, Ssym2

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j, it,jt, is,js, k,k1,k2, nd, ihm(ntrm),jhm(ntrm)
      Character(80) :: filename

      if(ncoef.le.0) Return

      write(filename,'(a,I3.3,I3.3)') 'idfc-res.',is,js
      open(3001,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'intc-res.',is,js
      open(3002,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'coef-res.',is,js
      open(3003,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'jhm.',is,js
!      open(3004,file=filename)
!      write(filename,'(a,I3.3,I3.3)') 'ijm.',is,js
!      open(3005,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'IPF-res.',is,js
      open(3006,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPF-res.',is,js
      open(3007,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPF-res.',is,js
      open(3008,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPF-res.',is,js
      open(3009,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'IPD-res.',is,js
      open(4006,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'JPD-res.',is,js
      open(4007,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'KPD-res.',is,js
      open(4008,file=filename)
      write(filename,'(a,I3.3,I3.3)') 'NPD-res.',is,js
      open(4009,file=filename)

!      write(3006,*) ndef,ldef
!      write(3006,*) IPF

!      write(3001,*) shape(idfc)
!      write(3001,*)
      write(3001,*) idfc(1:ncoef)

!      write(3002,*) shape(intc)
!      write(3002,*)
      write(3002,*) intc(1:ncoef)

!      write(3003,*) shape(coef)
!      write(3003,*)
      write(3003,*) coef(1:ntrm,1:ncoef)

      write(3006,*) ndef,ldef
      write(3006,*)
!      write(3006,*) shape(IPF) 
!      write(3006,*)
      write(3006,*) IPF(1:ndef)

!      write(3007,*) shape(JPF) 
!      write(3007,*)
      write(3007,*) JPF(1:ndef)

!      write(3008,*) shape(KPF) 
!      write(3008,*)
      write(3008,*) KPF(1:ndef)
      
!      write(3009,*) shape(NPF) 
!      write(3009,*)
      write(3009,*) NPF(1:ldef)

      write(4006,*) ndet,ldet
      write(4006,*)
!      write(3006,*) shape(IPF) 
!      write(3006,*)
      write(4006,*) IPD(1:ndet)

!      write(3007,*) shape(JPF) 
!      write(3007,*)
      write(4007,*) JPD(1:ndet)

!      write(3008,*) shape(KPF) 
!      write(3008,*)
      write(4008,*) KPD(1:ndet)

!      write(3009,*) shape(NPF) 
!      write(3009,*)
      write(4009,*) NPD(1:ndet)

! ... convert det.factors from ndef_list to common list:

      nd = 0; Do i=1,ncoef; nd=nd+idfc(i); End do
      if(nd.gt.0) then
       Call Ndet_Idet
       Do i=1,ncoef
        if(idfc(i).eq.0) Cycle
        j = idfc(i); idfc(i) = IPF(j)
        write(3006,*) i,j
       End do
      end if

!      write(3006,*) ndef,ldef,ncoef,nd
!      write(3006,*) 
!      write(3006,*) shape(IPF) 
!      write(3006,*)
!      write(3006,*) IPF

!      write(4006,*) ndet,ldet,ncoef,nd
!      write(4006,*)
!      write(3006,*) shape(IPF) 
!      write(3006,*)
!      write(4006,*) IPD

! ... define the term poiners:

      k = 0
      Do k1=1,kt1; it=IP_kt1(k1)
      Do k2=1,kt2; jt=IP_kt2(k2)
       if(is.eq.js.and.it.gt.jt) Cycle
       k = k + 1;  ihm(k) = it; jhm(k) = jt
      End do; End do

      if(k.ne.ntrm) Stop 'Add_res: ij <> ntrm'

! ... record the coef.s:

      if(.not.allocated(Cbuf)) &
      Allocate(Cbuf(mbuf),ibuf1(mbuf),ibuf2(mbuf),ibuf3(mbuf),ibuf4(mbuf))

      k = 0
      Do j = 1,ncoef
       Do i = 1,ntrm
        if(abs(coef(i,j)).lt.Eps_C) Cycle
        k = k + 1
        Cbuf(k)  = coef(i,j)
        ibuf1(k) = ihm(i)
        ibuf2(k) = jhm(i)
        ibuf3(k) = intc(j)
        ibuf4(k) = idfc(j)

        if(k.lt.mbuf) Cycle
        nbuf = k
        write(nu) nbuf
        write(nu) cbuf (1:nbuf)
        write(nu) ibuf1(1:nbuf)
        write(nu) ibuf2(1:nbuf)
        write(nu) ibuf3(1:nbuf)
        write(nu) ibuf4(1:nbuf)
!        write(3001,*) ibuf4(1:nbuf)
!        write(3002,*) ibuf3(1:nbuf)
!        write(3003,*) Cbuf(1:nbuf)
!        write(3004,*) ibuf2(1:nbuf)
!        write(3005,*) ibuf1(1:nbuf)
        k=0

       End do
      End do
      if(k.eq.0) Return

      nbuf = k
      write(nu) nbuf
      write(nu) cbuf (1:nbuf)
      write(nu) ibuf1(1:nbuf)
      write(nu) ibuf2(1:nbuf)
      write(nu) ibuf3(1:nbuf)
      write(nu) ibuf4(1:nbuf)
!      write(3001,*) ibuf4(1:nbuf)
!      write(3002,*) ibuf3(1:nbuf)
!      write(3003,*) Cbuf(1:nbuf)
!      write(3004,*) ibuf2(1:nbuf)
!      write(3005,*) ibuf1(1:nbuf)


      nc_new = nc_new + nbuf

      End Subroutine Add_res
