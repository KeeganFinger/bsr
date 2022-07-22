!======================================================================
      Subroutine Record_matrix_mpi(nu)
!======================================================================
!     collect and record of overlap matrix
!----------------------------------------------------------------------
      Use MPI
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,ip,jp,ich,jch, nn
      Real(8) :: S, v(ns), vv(ns,ns)
      Integer :: status(MPI_STATUS_SIZE)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... channel-channal blocks:

      nn = ns * ns
      print *, 'shape(hcc)', shape(hcc)
      Do ich = 1,nch
        Do jch = 1,ich
          i=icc(ich,jch)

          if(myid.ne.0.and.i.ne.0) then
            print *, 'sending',myid, ich, jch
            vv = hcc(1:ns,1:ns,i)
            Call MPI_SEND(vv,nn, MPI_DOUBLE_PRECISION, &
                         0, i, MPI_COMM_WORLD, ierr)
          end if

          if(myid.eq.0.and.i.ne.0) then
            print *, 'receiving',myid, ich, jch
            Call MPI_RECV(vv, nn, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                        i, MPI_COMM_WORLD, status, ierr)
            S = SUM(abs(vv))
            if(S.ne.0.d0) then
              write(nu) ich,jch
              write(nu) vv
            end if
          end if
          Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        End do
      End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

      Do ich = 1,nch
        Do ip = 1,npert
          i = icb(ich,ip)
          if(myid.ne.0.and.i.ne.0) then
            print *, 'sending',myid, ich, ip
            Call MPI_SEND(hcb(:,i),ns, MPI_DOUBLE_PRECISION, &
                    0, i, MPI_COMM_WORLD, ierr)
          endif
          if(myid.eq.0.and.i.ne.0) then
            print *, 'receiving',myid, ich, ip
            Call MPI_RECV(v, ns, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                    i, MPI_COMM_WORLD, status, ierr)
            S = SUM(abs(v))
            if(S.ne.0.d0) then
              write(nu) ip+nch,ich
              write(nu) v
            end if
          end if
          Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        End do
      End do

! ... perturter-perturber elements:

      Do ip = 1,npert
        Do jp = 1,ip
          i = ibb(ip,jp)
          if(myid.ne.0.and.i.ne.0) then
            print *, 'sending',myid, ip, jp
            Call MPI_SEND(hbb(i),1, MPI_DOUBLE_PRECISION, &
                        0, i, MPI_COMM_WORLD, ierr)
          endif
          if(myid.eq.0.and.i.ne.0) then
            print *, 'receiving',myid, ip, jp
            Call MPI_RECV(v, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                         i, MPI_COMM_WORLD, status, ierr)
            if(S.ne.0.d0) then
              write(nu) ip+nch,jp+nch
              write(nu) S
            end if
          end if
          Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        End do
      End do

      end if

! ... sign of the end

      if(myid.eq.0) write(nu) -1,-1

    End Subroutine Record_matrix_mpi
