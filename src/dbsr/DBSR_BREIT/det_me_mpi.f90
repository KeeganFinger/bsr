!======================================================================
      Subroutine DET_me
!======================================================================
!
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     m.e. between possible combinations of nj-orbitals
!
!----------------------------------------------------------------------

      USE nljm_orbitals;  USE conf_jj, ONLY: ne

      Implicit none
      Integer(4) :: i,i1,i2, j,j1,j2, k,k1,k2
      Integer(4), External :: Isort

!----------------------------------------------------------------------
!                              creat common list of orbital symmetries:
      ksym1=1; ksym2=1

      Nsym=0; k1=0;  k2=0; kz1=0; kz2=0

! ... exzaust the 1-st configuration:

      Do i = 1,ne 

       if(ksym1(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym1(i); Msym(Nsym)=Msym1(i); Jsym(Nsym)=Jsym1(i)
       k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=i; ksym1(i)=0

! ... check for the same orbitals rest the 1-st configuration:      
 
       Do j = i+1,ne
        if(ksym1(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym1(j)) Cycle
        if(Msym(Nsym).ne.Msym1(j)) Cycle
        if(Jsym(Nsym).ne.Jsym1(j)) Cycle
        k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=j; ksym1(j)=0
       End do        

!       Jdet=Isym1(1:ne);  kz1 = Isort(ne,Jdet) ???

! ... check for the same orbitals the 2-nd configuration:      

       IPsym2(Nsym)=k2
       Do j = 1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Jsym(Nsym).ne.Jsym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do        

      End do

      if(k1.ne.ne) Stop 'Det_me: k1 <> ne '

! ... exzaust the 2-st configuration:

      Do i = 1,ne 
       if(ksym2(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym2(i); Msym(Nsym)=Msym2(i); Jsym(Nsym)=Jsym2(i)
       k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=i; ksym2(i)=0
       IPsym1(Nsym)=k1

! ... check for the same orbitals rest of 2-st configuration:      
       
       Do j = i+1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Jsym(Nsym).ne.Jsym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do        

      End do

      if(k2.ne.ne) Stop 'Det_breit: k2 <> ne '

      Jdet=Isym1(1:ne);  kz1 = Isort(ne,Jdet) 
      Jdet=Isym2(1:ne);  kz2 = Isort(ne,Jdet)


!----------------------------------------------------------------------
!                              define the number of different orbitals:
      Ksym1(1)=ipsym1(1)
      Ksym2(1)=ipsym2(1)
      Do i = 2,NSYM
       Ksym1(i)=ipsym1(i)-ipsym1(i-1)
       Ksym2(i)=ipsym2(i)-ipsym2(i-1)
      End do

! ... how much different symmetries:

      k = 0
      Do i = 1,NSYM
       N1(i) = KSYM1(i)-KSYM2(i)
       N2(i) = KSYM2(i)-KSYM1(i)
       if(N1(i).gt.0) k = k + N1(i)
      End do

      if(k.gt.2) Return

!---------------------------------------------------------------------
!                                                         k = 2  case:
      Select case(k)

      Case(2)

       Do i=1,NSYM
        if(N1(i).le.0) Cycle; i1=i; N1(i)=N1(i)-1; Exit
       End do
       Do i=i1,NSYM
        if(N1(i).le.0) Cycle; i2=i; Exit
       End do
       Do i=1,NSYM
        if(N2(i).le.0) Cycle; j1=i; N2(i)=N2(i)-1; Exit
       End do
       Do i=j1,NSYM
        if(N2(i).le.0) Cycle; j2=i; Exit
       End do

       Call Zno_2ee(i1,i2,j1,j2)
  
!---------------------------------------------------------------------
!                                                         k = 1  case:
      Case(1)

       Do i=1,NSYM
        if(N1(i).le.0) Cycle; i1 = i; Exit
       End do
       Do i=1,NSYM
        if(N2(i).le.0) Cycle; j1 = i; Exit
       End do

       Do i = 1,Nsym
        if(Ksym1(i).eq.0) Cycle
        if(i.eq.i1.and.Ksym1(i).le.1) Cycle
        if(i.eq.j1.and.Ksym2(i).le.1) Cycle

        if(i.le.i1.and.i.le.j1)  then
          Call Zno_2ee(i,i1,i,j1)
        elseif(i.gt.i1.and.i.le.j1) then
          Call Zno_2ee(i1,i,i,j1)
        elseif(i.gt.i1.and.i.gt.j1) then
          Call Zno_2ee(i1,i,j1,i)
        elseif(i.le.i1.and.i.gt.j1) then
          Call Zno_2ee(i,i1,j1,i)
        end if

       End do

!---------------------------------------------------------------------
!                                                         k = 0  case:
      Case(0)

       Call ZNO_0ee
       Call ZNO_1ee

       Do i = 1,Nsym
        Do j = i,Nsym
         if(i.eq.j.and.Ksym1(i).le.1) Cycle;  Call Zno_2ee(i,j,i,j)
        End do
       End do

      End Select

      End Subroutine DET_me



!====================================================================
      Subroutine ZNO_0ee
!====================================================================
!
!     computes overlap integral between two determinants
!
!     Calls: Idet_fact, Incode_int, Iadd_zoef
!
!--------------------------------------------------------------------

      USE nljm_orbitals

      Implicit none
      Integer(4) :: idf,int
      Real(8) :: C
      Integer(4), External :: Idet_fact, Incode_int

      C = (-1)**(kz1+kz2)
      idf = Idet_fact (0,0,0,0)
      int = Incode_int (0,0,1,1,1,1)

      Call Iadd_zoef (C,int,idf)

      End Subroutine ZNO_0ee


!====================================================================
      Subroutine ZNO_1ee
!====================================================================
!
!    angular part of one-electron operator between two det.w.f
!
!    Calls: Idet_fact, Incode_int, Iadd_zoef.
!
!--------------------------------------------------------------------

      Use nljm_orbitals

      Implicit none

      Integer(4) :: i,j,i1,i2,k,k1,k2,is,idf,int
      Real(8) :: C
      Integer(4), External :: Idet_fact, Incode_int

      Do is=1,NSYM

       Do i=IPsym1(is-1)+1,IPsym1(is); k1=nnsym1(Isym1(i))
       Do j=IPsym1(is-1)+1,IPsym1(is); k2=nnsym2(Isym2(j))

       idf = Idet_fact(i,0,j,0)

       C=(-1)**(kz1+kz2+i+j)

       int = Incode_int (1,0,k1,k1,k2,k2)

       Call Iadd_zoef(C,int,idf)

       End do
       End do

      End do

      END Subroutine ZNO_1ee

!======================================================================
      SUBROUTINE ZNO_2ee (i1,i2,j1,j2)
!======================================================================
!
!     angular part of matrix elements between two det.w.f.
!     for two-electron operator
!
!     Calls: Check_boef, Idet_fact, Iadd_zoef.
!
!----------------------------------------------------------------------

      USE nljm_orbitals;  USE boef_list

      Implicit none
      Integer(4), Intent(in) :: i1,i2,j1,j2
      Integer(4) :: io1,io2,io3,io4
      Integer(4) :: i,k, k1,k2,k3,k4, int,idf,kz

      Integer(4), External :: Idet_fact, Incode_int

!----------------------------------------------------------------------

      if(mod(Lsym(i1)+Lsym(i2)+Lsym(j1)+Lsym(j2),2).ne.0) Return

      if(Msym(i1)+Msym(i2).ne.Msym(j1)+Msym(j2)) Return

!----------------------------------------------------------------------

      Call Check_boef(Lsym(i1),Jsym(i1),Msym(i1), & 
                      Lsym(i2),Jsym(i2),Msym(i2), &
                      Lsym(j1),Jsym(j1),Msym(j1), &
                      Lsym(j2),Jsym(j2),Msym(j2))

!----------------------------------------------------------------------
!
      Do k1 = IPsym1(i1-1)+1,IPsym1(i1); io1=nnsym1(Isym1(k1))
      Do k2 = IPsym1(i2-1)+1,IPsym1(i2); io2=nnsym1(Isym1(k2))  
       if(k2.le.k1) Cycle 
                       
      Do k3 = IPsym2(j1-1)+1,IPsym2(j1); io3=nnsym2(Isym2(k3))
      Do k4 = IPsym2(j2-1)+1,IPsym2(j2); io4=nnsym2(Isym2(k4))  
       if(k4.le.k3) Cycle 

       idf = Idet_fact(k1,k2,k3,k4)

       kz = (-1)**(kz1+kz2+k1+k2+k3+k4)

       Do i = ncblk(kblk-1)+1,ncblk(kblk)
        int = IB_int(i)
        if(int.gt.0) then
         k = int-1
         int = Incode_int(2,k,io1,io2,io3,io4)
         Call Iadd_zoef(Boef(i)*kz,int,idf)
        else
         k = -int-1
         int = Incode_int(2,k,io1,io2,io4,io3)
         Call Iadd_zoef(Boef(i)*kz,int,idf)
        end if        
       End do

      End do;  End do;  End do;  End do

      END SUBROUTINE ZNO_2ee


!======================================================================
      Integer(4) Function Idet_fact (i1,i2,i3,i4)
!======================================================================
!
!     determines the overlap factor and its position in NDEF list
!     for matrix element between two determinant wave functions
!     located in module 'nljm_orbitals' 
!
!     (j1,j2) and (j3,j4) - active electrons in the first and second 
!                           determinants 
!
!     Calls:  Iadd_ndet, Iadd_ndef, ISORT
!
!----------------------------------------------------------------------

      USE param_jj, ONLY: ibd,ibf
      USE nljm_orbitals

      Implicit none
      Integer(4), Intent(in) :: i1,i2,i3,i4
      Integer(4) :: i, k,k1,k2, is,id,kd
      Integer(4), External :: Iadd_ndet, Iadd_ndef, ISORT

      kd = 0
      Do is = 1,NSYM

       k1 = 0
       Do i = IPsym1(is-1)+1,IPsym1(is)
        if(i.eq.i1.or.i.eq.i2) Cycle
        k1 = k1 + 1;  N1(k1) = nnsym1(Isym1(i))
       End do

       k2 = 0
       Do i = IPsym2(is-1)+1,IPsym2(is)
        if(i.eq.i3.or.i.eq.i4) Cycle
        k2 = k2 + 1;  N2(k2) = nnsym2(Isym2(i))
       End do

       if(k1.ne.k2) Stop 'Idet_fact: k1 <> k2'

       if(k1.eq.0) Cycle

       NP(1:k1) = N1(1:k1)*ibd + N2(1:k1);  id = Iadd_ndet(k1,NP)

       Do k=1,kd
        if(N3(k).ne.id) Cycle;  N4(k)=N4(k)+1; id=0; Exit
       End do
       if(id.eq.0) Cycle

       kd=kd+1; N3(kd)=id; N4(kd)=1

      End do

      Idet_fact=0; if(kd.eq.0) Return

      NP(1:kd)=N3(1:kd)*ibf + N4(1:kd); k=ISORT(kd,NP)

      Idet_fact = Iadd_ndef(kd,NP)

      End Function Idet_fact


!======================================================================
      Integer(4) Function Incode_int(met,k,I1,I2,I3,I4)
!======================================================================
!
!     incode integral
!
!----------------------------------------------------------------------
 
      Use param_jj

      Implicit none
      Integer(4), intent(in) :: met,k,I1,I2,I3,I4       

      if(max(i1,i2,i3,i4).ge.ib5) Stop ' INT_pack: i > pack-base '

      Incode_INT = ((((i1*ib5+i2)*ib5+i3)*ib5+i4)*ib10+k)*ib2+met

      End Function Incode_INT


!======================================================================
      Subroutine Decode_int (met,k,I1,I2,I3,I4,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use param_jj

      Implicit none
      Integer(4) :: int, met,k,I1,I2,I3,I4, ii

      ii = int
      met = mod(ii,ib2);  ii = ii/ib2
      k   = mod(ii,ib10); ii = ii/ib10
      I4  = mod(ii,ib5);  ii = ii/ib5
      I3  = mod(ii,ib5);  ii = ii/ib5
      I2  = mod(ii,ib5);  I1 = ii/ib5

      End Subroutine Decode_INT
