!======================================================================
      Module radial_overlaps
!======================================================================
!     contains the desription of one-electron radial overlaps 
!----------------------------------------------------------------------
      Implicit none
   
      Integer :: mobs = 0      !  max.number of  overlaps
!      Integer :: nobs = 0      !  curent number of overlaps
      Integer :: kobs = 2**10  !  initial suggestion for mobs

!     value of the bound-bound one-electron overlap:

!      Real(8), allocatable :: Cobs(:)

!     pointer to orbitals (iobs > = jobs):

      Integer, allocatable :: iobs(:),jobs(:)

      End Module radial_overlaps

 
!======================================================================
      Subroutine Alloc_radial_overlaps(m)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use radial_overlaps
      Use orb_overlaps

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
      Use radial_overlaps
      Use orb_overlaps

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
      Use radial_overlaps
      Use orb_overlaps

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
