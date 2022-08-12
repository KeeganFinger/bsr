!======================================================================
      Subroutine Alloc_orb_overlaps(nbf,lbs,iech,ncore)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use orb_overlaps

      Implicit none
      Integer, intent(in) :: nbf,lbs(nbf),iech(nbf),ncore
      Real(8), external :: QUADR
      Integer :: i,j,l, il,jl,is,js, ip 

      if(allocated(ipl)) Deallocate(ipl,jpl,lorb,chan,ip_l,jp_l,ip_ovl,Cobs)
      norb = 0; mem_orb_overlaps = 0
      if(nbf.le.0) Return

! ... sort the orbitals according l-values:
  
      norb = nbf
      Allocate(lorb(norb)); lorb(1:norb) = lbs(1:norb)
      Allocate(chan(norb)); chan(1:norb) = iech(1:norb)
      Allocate(ipl(norb),jpl(norb))
      Call SORTI(norb,lorb,ipl)

      Do i=1,norb; j = ipl(i); jpl(j) = i; End do

! ... find last entry for each l:

      max_l  = maxval(lorb)
      Allocate(ip_l(0:max_l),jp_l(0:max_l));  jp_l=0
      Do j = 1, norb; i = ipl(j)
       l=lorb(i); jp_l(l) = j
      End do

! ... ip_l is  "base" for given l:     

      Do l=1,max_l;  if(jp_l(l).eq.0) jp_l(l)=jp_l(l-1); End do
      Do l=max_l,1,-1;  ip_l(l) = jp_l(l-1); End do;  ip_l(0)=0

! ... allocate the overlaps:

      Allocate(ip_ovl(0:max_l)); ip_ovl=0
      i = jp_l(0)-ip_l(0); nobs = i*(i+1)/2  
      Do l = 1,max_l
       ip_ovl(l) = nobs
       i = jp_l(l)-ip_l(l)
       nobs = nobs + i*(i+1)/2  
      End do
      Allocate(Cobs(nobs)); Cobs = 0.d0
      mem_orb_overlaps = 2*nobs + 4*norb + 3*max_l

! ... radial overlaps:

!     il,jl - indexes in the ordered list

      Do il = 1,norb;  i = ipl(il);  l = lorb(i)
      Do jl = 1,il;    j = ipl(jl);  if(l.ne.lorb(j)) Cycle              

       ip = ip_l(l)
       is = max(il,jl)-ip; js = min(il,jl)-ip
       ip = ip_ovl(l) + is*(is-1)/2 + js

       if(chan(i).eq.0.and.chan(j).eq.0) then
        Cobs(ip) = QUADR(i,j,0)
       elseif(chan(i).ne.0.and.chan(j).eq.0) then
        Cobs(ip) = i*ibo+j        
        if(j.le.ncore) Cobs(ip) = 0.d0
       elseif(chan(i).eq.0.and.chan(j).ne.0) then
        Cobs(ip) = j*ibo+i        
        if(i.le.ncore) Cobs(ip) = 0.d0
       else
        Cobs(ip) = i*ibo+j        
       end if

      End do
      End do

      Call Check_orb_overlaps

      End Subroutine Alloc_orb_overlaps


