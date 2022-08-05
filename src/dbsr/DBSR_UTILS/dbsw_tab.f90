!======================================================================
!     UTILITY       D B S W _ T A B
!
!               C O P Y R I G H T -- 2008
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     connverts the B-spline radial orbital wbs-files into tab-files
!     suitable for grafic display
!----------------------------------------------------------------------
!
!     INPUT FILE:    name.bsw
!                    
!     OUTPUT FILES:  name.bsw.nl  for each nl
!
!----------------------------------------------------------------------
!     ARGUMENTS:     name.bsw  
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq
      Use zconst, only: c_au

      Implicit real(8) (A-H,O-Z)
      Character(1) :: ans
      Character(5) :: EL
      Character(40) :: AF,BF
      Real(8), allocatable ::  R(:),P(:),Q(:)
      Real(8) :: nuc_charge

! ... input data: 
         
      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,AF)

      if(iarg.lt.1.or.AF.eq.'?') then
       write(*,*)
       write(*,*) 'dbsw_tab connverts the B-spline radial orbital wbs-files into'
       write(*,*) 'text files suitable for grafic display'
       write(*,*)
       write(*,*) 'Call as:   dbsw_tab  name.bsw'
       write(*,*)
       write(*,*) 'file will be created for each orbital'
       Stop ' '
      end if

      if(iarg.eq.3) then
        Call read_rarg('nuc_charge',nuc_charge)
        Call read_rarg('awt',awt)
      endif

! ... set up B-splines:
 
      Call Check_file('knot.dat')
      Call def_grid('knot.dat',AF,nuc_charge,awt)
      Call alloc_DBS_gauss
!      Call alloc_DBS_galerkin

! ... radial w.f.:

      nuw=1
      Open(nuw,file=AF,status='OLD',form='UNFORMATTED')
      Call Read_pqbs(nuw)
      Close(nuw)
      iaf = LEN_TRIM(AF)

! ... sets up grid points and initializes the values of the spline: 

      NR = nv*ks+2; Allocate(R(NR),P(NR),Q(NR))
      ii=1; R(1)=0.d0
      Do i=1,nv; Do j=1,ks; ii=ii+1; R(ii) = gr(i,j); End do; End do
      ii=ii+1; R(ii) = t(ns+1)

! ... Cycle over nl in input:

      BF=AF; iaf=iaf+1; BF(iaf:iaf)='.'

      Do i=1,nbf

       ip=i+i-1; jp=ip+1
       P=0.d0; Call Bvalue_bm(ksp,pq(1,1,i),P,pbsp)
       Q=0.d0; Call Bvalue_bm(ksq,pq(1,2,i),Q,qbsp)

       k=iaf; EL='     '; EL=ebs(i)
       Do j=1,5
        if(EL(j:j).eq.' ') Cycle
        k=k+1; BF(k:k)=EL(j:j)
       End do
       iout=2; Open(iout,file=BF)
       write(iout,'(3(5x,a,10x))') 'R','P','Q'

       S=1.d0
!       S = c_au * 2.d0 * t(ns+1) / kbs(i)

       Do j=1,nr
        write(iout,'(3D16.8)') R(j),P(j),Q(j)*S
       End do

      End do

      
      END   !  program dbsr_tab
          

