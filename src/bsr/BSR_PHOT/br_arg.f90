!======================================================================
      Subroutine br_arg
!======================================================================

      USE MPI
      Use bsr_phot

      Implicit none

      Call MPI_BCAST(klsp,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iauto, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mfgi,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwt,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ikm,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ibug,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(AC,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(DR,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(AWT,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(e_exp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(athreshold,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(bthreshold,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(me,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) Allocate(EK(me),IEK(me))

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(iek,me,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ek,me,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg
