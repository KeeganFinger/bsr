!======================================================================
      Module overlaps
!======================================================================
!     contains the desription of one-electron radial overlaps
!----------------------------------------------------------------------
      Implicit none

      Integer :: mobs = 0      !  max.number of  overlaps
      Integer :: nobs = 0      !  curent number of overlaps
      Integer :: kobs = 2**10  !  initial suggestion for mobs
      Integer :: norb = 0       ! number of orbitals

! ... pointer to orbitals in l-order, (1:norb) arrays:

      Integer, allocatable :: lorb(:)   ! list of l-values for all orbitals
      Integer, allocatable :: ipl(:)    ! l-values ordering pointer
      Integer, allocatable :: jpl(:)    ! jpl(ipl(i)) = i

! ... base pointer to orbitals for given l, (0:max_l) arrays

      Integer :: max_l = 0
      Integer, allocatable :: ip_l(:)   ! first -1 entry (base) for given l in ordered list, ip_l(l)=jp_l(l-1)
      Integer, allocatable :: jp_l(:)   ! last entry for given l in ordered list

! ... base pointer to overlaps for given l in Cobs array,
! ... where all overlaps are presented half square-matrixes

      Integer, allocatable :: ip_ovl(:)

! ... value of the bound-bound one-electron overlap
! ... or orthogonal conditions for continuum orbitals

      Real(8), allocatable :: Cobs(:)

! ... memory in words (4 bytes):

      Integer :: mem_orb_overlaps = 0

      Integer :: ibo = 2**15

      Integer, allocatable :: iobs(:),jobs(:)

      End Module overlaps


!======================================================================
      Subroutine Alloc_radial_overlaps(m)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use overlaps

      Implicit none
      Integer, intent(in) :: m
      Real(8), allocatable :: rarray(:)
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(allocated(Cobs)) Deallocate(Cobs,iobs,jobs)
       mobs = kobs; nobs = 0
       Allocate(Cobs(mobs),iobs(mobs),jobs(mobs))
      elseif(m.eq.0) then
       if(allocated(Cobs)) Deallocate(Cobs,iobs,jobs)
       mobs = 0; nobs = 0
      elseif(m.gt.mobs) then
       if(nobs.le.0) then
        if(Allocated(Cobs)) Deallocate(Cobs,iobs,jobs)
        mobs = m; nobs = 0
        Allocate(Cobs(mobs),iobs(mobs),jobs(mobs))
       else
        Allocate(rarray(nobs)); rarray(1:nobs)=Cobs(1:nobs)
        Deallocate(Cobs); Allocate(Cobs(m))
        Cobs(1:nobs)=rarray(1:nobs); Deallocate(rarray)
        Allocate(iarray(nobs)); iarray(1:nobs)=iobs(1:nobs)
        Deallocate(iobs); Allocate(iobs(m))
        iobs(1:nobs)=iarray(1:nobs)
        iarray(1:nobs)=jobs(1:nobs)
        Deallocate(jobs); Allocate(jobs(m))
        jobs(1:nobs)=iarray(1:nobs); Deallocate(iarray)
        mobs = m
       end if
      end if

      End Subroutine Alloc_radial_overlaps


!======================================================================
      Subroutine Iadd_obs(io,jo,S)
!======================================================================
!     add new entry( or replace) in the ordered list of overlaps
!----------------------------------------------------------------------
      Use overlaps

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S
      Integer :: i,j, k,l,m

      if(mobs.eq.0) Call Alloc_radial_overlaps(kobs)

      i = max(io,jo)
      j = min(io,jo)

! ... search position (m) for given overlap

      k=1; l=nobs
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (i.lt.iobs(m)) then;     l = m - 1
      elseif(i.gt.iobs(m)) then;     k = m + 1
      else
       if    (j.lt.jobs(m)) then;    l = m - 1
       elseif(j.gt.jobs(m)) then;    k = m + 1
       else
        Cobs(m)=S;  Return
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:

      Do l=nobs,k,-1; m = l + 1
       Cobs(m)=Cobs(l); iobs(m)=iobs(l); jobs(m)=jobs(l)
      End do

! ... add new overlap:

      Cobs(k)=S; iobs(k)=i; jobs(k)=j; nobs=nobs+1
      if(nobs.eq.mobs) Call alloc_radial_overlaps(mobs+kobs)

      End Subroutine Iadd_obs


!======================================================================
      Subroutine Idel_obs(io)
!======================================================================
!     add new entry( or replace) in the ordered list of overlaps
!----------------------------------------------------------------------
      Use overlaps

      Implicit none
      Integer, intent(in) :: io
      Integer :: i,k

      k = 0
      Do i=1,nobs
       if(io.eq.iobs(i).or.io.eq.jobs(i)) Cycle
       k = k+1; if(k.eq.i) Cycle
       iobs(k)=iobs(i); jobs(k)=jobs(i)
      End do
      nobs = k

      End Subroutine Idel_obs

! !======================================================================
!       Real(8) Function OBS(io,jo)
! !======================================================================
! !     find overlaps <io|jo>
! !----------------------------------------------------------------------
!       Use overlaps
!
!       Implicit none
!       Integer, intent(in) :: io,jo
!       Integer :: i,j, k,l,m
!
!       i = max(io,jo)
!       j = min(io,jo)
!
! ! ... search position (m) for given overlap
!
!       k=1; l=nobs
!     1 if(k.gt.l) go to 2
!       m=(k+l)/2
!       if    (i.lt.iobs(m)) then;     l = m - 1
!       elseif(i.gt.iobs(m)) then;     k = m + 1
!       else
!        if    (j.lt.jobs(m)) then;    l = m - 1
!        elseif(j.gt.jobs(m)) then;    k = m + 1
!        else
!         OBS=Cobs(m);  Return
!        end if
!       end if
!       go to 1
!     2 OBS = 0.d0
!
!       End Function OBS

!======================================================================
      Real(8) Function OBS(io,jo)
!======================================================================
!     find overlaps <io|jo>
!----------------------------------------------------------------------
      Use overlaps

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j,l, ip, is,js

! ... check the arguments:

      obs = 0.d0; if(lorb(io).ne.lorb(jo)) Return;  l=lorb(io)
      i = jpl(io);  j = jpl(jo)
      ip = ip_l(l)
      is = max(i,j)-ip; js = min(i,j)-ip
      ip = ip_ovl(l) + is*(is-1)/2 + js

      obs = Cobs(ip)

      End Function OBS


!======================================================================
      Real(8) Function VDET (kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant for given orbitals
!
!     Calls:  DET
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: kd,N1(kd),N2(kd)
      Integer :: i,j
      Real(8) :: ADET(kd*kd)
      Real(8), external :: DET, OBS

      if(kd.eq.0) then
       VDET = 1.d0
      elseif(kd.eq.1) then
       VDET = OBS(N1(1),N2(1))
      elseif(kd.eq.2) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))
      elseif(kd.eq.3) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2))*OBS(N1(3),N2(3)) +  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(3))*OBS(N1(3),N2(1)) +  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(1))*OBS(N1(3),N2(2)) -  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(2))*OBS(N1(3),N2(1)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))*OBS(N1(3),N2(3)) -  &
              OBS(N1(1),N2(1))*OBS(N1(2),N2(3))*OBS(N1(3),N2(2))
      else
       Do i=1,kd;  Do j=1,kd
         adet((i-1)*kd+j)=OBS(N1(i),N2(j))
       End do; End do
       VDET = DET(kd,adet)
      end if

      End Function VDET


!======================================================================
      Subroutine Check_orb_overlaps
!======================================================================
!     the imposed orthogonal conditions given as < n1k1| n2k2>=0
!----------------------------------------------------------------------
      Use overlaps
      Use internal_file

      Implicit none
      Integer :: i,j,l, i1,n1,l1,k1, i2,n2,l2,k2, is,js,ip, ii
      Integer, external :: Ifind_bsorb
      Character(13) :: Aort

      Do ii = 1,nlines

      Aort = trim(aline(ii))
      if(Aort(1:1).ne.'<') Cycle
      Call EL4_nlk(Aort(2: 5),n1,l1,k1)
      i1 = Ifind_bsorb(n1,l1,k1)
      if(i1.eq.0) Cycle
      Call EL4_nlk(Aort(7:10),n2,l2,k2)
      i2 = Ifind_bsorb(n2,l2,k2)
      if(i2.eq.0) Cycle

      if(l1.ne.l2) Cycle; l = l1
      read(Aort(13:13),'(i1)') j
      if(j.ne.0) Cycle

      i = jpl(i1);  j = jpl(i2)
      ip = ip_l(l)
      is = max(i,j)-ip; js = min(i,j)-ip
      ip = ip_ovl(l) + is*(is-1)/2 + js

      Cobs(ip) = 0.d0

      End do

      End Subroutine Check_orb_overlaps
