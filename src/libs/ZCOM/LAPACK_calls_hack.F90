!----------------------------------------------------------------------
!  interface routines to call some frequently-used LAPACK routines
!
!  LAP_DGESV     -    A x = B,      A,B - real
!  LAP_ZGESV     -    A x = B,      A,B - complex
!  LAP_DSYEV     -    A x = E x,    A - real symmetric
!  LAP_DSYEVX    -    A x = E x,    for some eigenvalues only
!  LAP_ZHEEV     -    A x = E x,    A - complex Hermitian
!  LAP_DSYGV     -    A x = E C x,  A,C - real symmetric
!  LAP_DSYGVX    -    A x = E C x,  for some eigenvalues only
!  LAP_INV       -    A^-1,         A - real
!  LAP_INVS      -    A^-1,         A - real symmetric
!----------------------------------------------------------------------

!======================================================================
      Subroutine LAP_DPOTRF(UPLO,n,A,lda,info)
!======================================================================
!     Call LAPACK procedure DPORTF to compute the Cholesky
!     factorization of a real symmetric positive definite matrix
!     A(lda,n)
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!----------------------------------------------------------------------
!@CUF USE CUSOLVERDN
      Implicit none
      Character(1), intent(in) :: UPLO
      Integer, intent(in) :: n,lda
      Integer, intent(out) :: info
      Real(8) :: A(lda,n)
!@CUF ATTRIBUTES(MANAGED) :: A

#ifdef _CUDA
      type(cusolverDnHandle) :: handle
      Real(8), allocatable, managed :: work(:)
      Integer :: istat, ijobz, iuplo, lwork
      Integer, device :: devInfo(1)
      Integer :: cusolverConvertToFillMode

      istat = cusolverDnCreate(handle)
      iuplo = cusolverConvertToFillMode(UPLO)
      istat = cusolverDnDpotrf_buffersize(handle, iuplo, &
                     n, A, lda, lwork)
      if (istat.ne.0) write(*,*) &
        'LAP_DSYEV: DSYEV(lapack) cusolver gives ISTAT = ',istat

      allocate(work(lwork))

      istat = cusolverDnDpotrf(handle, iuplo, n, A, lda, &
                    work, lwork, devInfo(1))
      if (istat.ne.0) write(*,*) &
        'LAP_DSYEV: DSYEV(lapack) cusolver gives ISTAT = ',istat

      INFO = devInfo(1)
      deallocate(work)

#else
      Call DPOTRF(UPLO,n,A,lda,INFO)

#endif
      if(INFO.ne.0) write(*,*) 'LAP_DSYEV: DSYEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_DPOTRF

!======================================================================
      Subroutine LAP_DGESV(m,n,k,A,B,info)
!======================================================================
!     Call LAPACK procedure DGESV to solve the system of algebraic
!     equations A x = B, results in B.
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m,k
      Real(8) :: A(m,*),B(m,*)
      Integer, intent(out), optional :: info
      Integer :: IPIV(m)

      Call DGESV( n,k, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' DGESV(lapack) give INFO =',INFO

      End Subroutine LAP_DGESV


!======================================================================
      Subroutine LAP_ZGESV(n,m,A,B,info)
!======================================================================
!     Call LAPACK procedure ZGESV to solve the system of complex
!     algebraic  equations A x = B, results in B.
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      COMPLEX(8) :: A(m,m),B(m,m)
      Integer, intent(out), optional :: info
      Integer :: IPIV(m)

      Call ZGESV( m,n, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' ZGESV(lapack) give INFO =',INFO

      End Subroutine LAP_ZGESV


!======================================================================
      Subroutine LAP_DSYEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for Real symmetric matrix A(n,n)
!     job = 'V' or 'N' - compute or not the eigenvectors
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
!@CUF USE CUSOLVERDN
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out) :: info
      Real(8) :: A(m,m)
      Real(8) :: eval(m)
!@CUF ATTRIBUTES(MANAGED) :: A, eval

#ifdef _CUDA
      type(cusolverDnHandle) :: handle
      Real(8), allocatable, managed :: work(:)
      Integer :: istat, ijobz, iuplo, lwork
      Integer, device :: devInfo(1)
      Integer :: cusolverConvertToFillMode
      Integer :: cusolverConvertToEigMode

      istat = cusolverDnCreate(handle)
      ijobz = cusolverConvertToEigMode(job)
      iuplo = cusolverConvertToFillMode(UPLO)
      istat = cusolverDnDsyevd_buffersize(handle, ijobz, iuplo, &
                     N, A, N, eval, lwork)
      if (istat.ne.0) write(*,*) &
        'LAP_DSYEV: DSYEV(lapack) cusolver gives ISTAT = ',istat

      allocate(work(lwork))

      istat = cusolverDnDsyevd(handle, ijobz, iuplo, N, A, N,  &
                     eval, work, lwork, devInfo(1))
      if (istat.ne.0) write(*,*) &
        'LAP_DSYEV: DSYEV(lapack) cusolver gives ISTAT = ',istat

      INFO = devInfo(1)
      deallocate(work)

#else
      Real(8) :: work(3*n-1)
      Integer :: lwork

      lwork = 3*n-1

      Call DSYEV(job,UPLO,n,A,m,eval,WORK,LWORK,INFO )

#endif
      if(INFO.ne.0) write(*,*) 'LAP_DSYEV: DSYEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_DSYEV


!======================================================================
      Subroutine LAP_DSYEVX(job,UPLO,n,m,A,eval,k,info)
!======================================================================
!     Call LAPACK procedure DSYEVX to obtain the selected eigenvalues
!     (eval) and eigenvectors (A) for problem A S = E S
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m,k
      Integer, intent(out), optional :: info
      Real(8) :: A(m,*)
      Real(8) :: eval(*)

      Real(8), allocatable :: work(:), Z(:,:)
      Integer, allocatable :: iwork(:),ifail(:)
      Integer :: kk, lwork
      Real(8) :: ABSTOL, VL =0.d0, VU = 0.d0
      Real(8), external :: DLAMCH

      ABSTOL = 4*DLAMCH('S')
      lwork = 10*n
      Allocate(work(lwork),iwork(5*n),ifail(k),Z(n,k))

      Call DSYEVX(job,'I',UPLO,n,A,m,VL,VU,1,k,ABSTOL,kk,eval, &
                  Z,n, WORK, LWORK, IWORK, IFAIL, INFO)

      if(kk.ne.k) then
       write(*,*) ' DSYEVX(lapack) provides',kk,' eigenvectors'
       write(*,*) ' when we ordered', k,'  ones'
       Stop ' Stop in LAP_DSYEVX'
      end if

      if(INFO.ne.0) write(*,*) 'DSYEVX(lapack) gives INFO = ',INFO

      A(1:n,1:k) = Z(1:n,1:k)

      Deallocate(work,iwork,ifail,Z)

      End Subroutine LAP_DSYEVX


!======================================================================
      Subroutine LAP_DSYGV(job,UPLO,m,A,lda,C,ldc,eval,info)
!======================================================================
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for generalized problem A S = E C S
!     job = 'N' - compute eigenvalues only
!           'V' - compute eigenvalues and eigenvectors
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
#ifdef _CUDA
      use cusolverdn
#endif
!!@CUF USE CUSOLVER
            Implicit none
            Character(1), intent(in) :: job,UPLO
            Integer, intent(in) :: m,lda,ldc
            Integer, intent(out), optional :: info
            Real(8) :: A(lda,*),C(ldc,*)
            Real(8) :: eval(m)
!!@CUF ATTRIBUTES(MANAGED) :: A, C, eval
#ifdef _CUDA
      ATTRIBUTES(MANAGED) :: A, C, eval
#endif

      Integer :: lwork

#ifdef _CUDA
      type(cusolverDnHandle) :: handle
      Integer :: cusolverConvertToEigMode, cusolverConvertToFillMode
      Real(8), allocatable, managed :: work(:)
      Integer :: istat, ijobz, iuplo
      Integer, device :: devInfo(1)

      istat = cusolverDnCreate(handle)
      ijobz = cusolverConvertToEigMode(job)
      iuplo = cusolverConvertToFillMode(UPLO)

      istat = cusolverDnDsygvd_buffersize(handle, 1, ijobz, iuplo, &
                        m, A, lda, C, ldc, eval, lwork)
      if (istat.ne.0) write(*,*) &
            'LAP_DSYGV: DSYGV(lapack) cusolver gives ISTAT = ',istat

      allocate(work(lwork))

      istat = cusolverDnDsygvd(handle, 1, ijobz, iuplo, m, A, lda, C, &
                ldc, eval, work, lwork, devInfo(1))
      if (istat.ne.0) write(*,*) &
            'LAP_DSYGV: DSYGV(lapack) cusolver gives ISTAT = ',istat

      INFO = devInfo(1)
      deallocate(work)
#else
      Real(8) :: work(3*m)

      lwork = 3*m

      Call DSYGV(1,job,UPLO,m,A,lda,C,ldc,eval,WORK,LWORK,INFO)

      if(INFO.ne.0) write(*,*) ' DSYGV(lapack) gives INFO = ',INFO
#endif


      End Subroutine LAP_DSYGV


!======================================================================
      Subroutine LAP_DSYGVX(job,UPLO,n,m,A,C,eval,k,info)
!======================================================================
!     Call LAPACK procedure DSYGVX to obtain first k eigenvalues
!     (eval) and eigenvectors (A) for generalized problem  A S = E C S
!
!     job = 'N' - compute eigenvalues only
!           'V' - compute eigenvalues and eigenvectors
!
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m,k
      Integer, intent(out), optional :: info
      Real(8) :: A(m,*),C(m,*)
      Real(8) :: eval(*)

      Real(8), allocatable :: work(:), Z(:,:)
      Integer, allocatable :: iwork(:),ifail(:)

      Integer :: kk, lwork
      Real(8) :: ABSTOL, VL =0.d0, VU = 0.d0
      Real(8), external :: DLAMCH

      ABSTOL = 4*DLAMCH('S')
      lwork = 10*n
      Allocate(work(lwork),iwork(5*n),ifail(k),Z(n,k))

      Call DSYGVX(1,job,'I',UPLO,n,A,m,C,m,VL,VU,1,k,ABSTOL,kk,eval, &
                  Z,n, WORK, LWORK, IWORK, IFAIL, INFO)

      if(kk.ne.k) then
       write(*,*) ' DSYGVX(lapack) provides',kk,' eigenvectors'
       write(*,*) ' when we ordered', k
      end if
      if(INFO.ne.0) write(*,*) 'DSYGVX(lapack) gives INFO = ',INFO

      A(1:n,1:k) = Z(1:n,1:k)

      Deallocate(work,iwork,ifail,Z)

      End Subroutine LAP_DSYGVX


!======================================================================
      Subroutine LAP_INV(n,m,A,info)
!======================================================================
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for matrix A
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      Real(8) :: A(m,m)
      Integer, intent(out), optional :: info
      Real(8) :: work(3*n)
      Integer :: ipiv(n), lwork,ierr

      Call DGETRF(n,n,A,m,IPIV,info)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO

      lwork = 3*n
      Call DGETRI(n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      End Subroutine LAP_INV


!======================================================================
      Subroutine LAP_INVS(UPLO,n,m,A,info)
!======================================================================
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for symmetric matrix
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      Character(1), intent(in) :: UPLO
      Real(8) :: A(m,m)
      Integer, intent(out), optional :: info
      Real(8) :: work(3*n)
      Integer :: ipiv(n), lwork,ierr

      lwork = 3*n

      Call DSYTRF(UPLO,n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO
      if(info.ne.0) Return

      Call DSYTRI(UPLO,n,A,m,IPIV,WORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      End Subroutine LAP_INVS


!======================================================================
      Subroutine LAP_ZHEEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure ZHEEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for complex Hermitian matrix A(n,n)
!     job = 'V' or 'N' - compute or not the eigenvectors
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out), optional :: info
      COMPLEX(8) :: A(m,m)
      Real(8) :: eval(m)
      Real(8) :: rwork(3*n)
      COMPLEX(8) :: work(2*n)
      Integer :: lwork,lrwork

      lwork = 2*n; lrwork= 3*n

      Call ZHEEV(job,UPLO,n,A,m,eval,WORK,lwork,RWORK,INFO)

      if(INFO.ne.0) write(*,*) ' ZHEEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_ZHEEV

      SUBROUTINE LAP_DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!@CUF USE CUBLAS
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!@CUF ATTRIBUTES(MANAGED) :: A, B
!     ..
      Call dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
      end subroutine
!
!  =====================================================================
!
      SUBROUTINE LAP_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!@CUF USE CUBLAS
!@CUF USE CUDAFOR
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!@CUF ATTRIBUTES(MANAGED) :: A, B, C
!     ..
      Call dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
!@CUF istat = cudaDeviceSynchronize()
      end subroutine

!  =====================================================================

      SUBROUTINE LAP_DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!@CUF USE CUBLAS
!@CUF USE CUDAFOR
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!@CUF ATTRIBUTES(MANAGED) :: A, B
!
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER, LOWER
      INTEGER            K, KB, NB

!@CUF ISTAT = cudaMemAdvise(B, LDB*N, cudaMemAdviseSetReadMostly, 0)

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LOWER = LSAME( UPLO, 'L' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT. UPPER .AND. .NOT. LOWER ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGST', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
      NB = 64
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
         CALL DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
!
!        Use blocked code
!
        IF( ITYPE.EQ.1 ) THEN
          IF( UPPER ) THEN
!
!           Compute inv(U**T)*A*inv(U)
!
            DO 10 K = 1, N, NB
              KB = MIN( N-K+1, NB )
!
!             Update the upper triangle of A(k:n,k:n)
!
              CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                           B( K, K ), LDB, INFO )
              IF( K+KB.LE.N ) THEN
                 CALL DTRSM( 'Left', UPLO, 'Transpose', 'Non-unit', &
                             KB, N-K-KB+1, ONE, B( K, K ), LDB, &
                             A( K, K+KB ), LDA )
                 CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                             A( K, K ), LDA, B( K, K+KB ), LDB, ONE, &
                             A( K, K+KB ), LDA )
                 CALL DSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE, &
                              A( K, K+KB ), LDA, B( K, K+KB ), LDB, &
                              ONE, A( K+KB, K+KB ), LDA )
                 CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                             A( K, K ), LDA, B( K, K+KB ), LDB, ONE, &
                             A( K, K+KB ), LDA )
                 CALL DTRSM( 'Right', UPLO, 'No transpose', &
                             'Non-unit', KB, N-K-KB+1, ONE, &
                             B( K+KB, K+KB ), LDB, A( K, K+KB ), &
                             LDA )
              END IF
!@CUF istat = cudaDeviceSynchronize()
   10       CONTINUE
          ELSE
!
!           Compute inv(L)*A*inv(L**T)
!
            DO 20 K = 1, N, NB
              KB = MIN( N-K+1, NB )
!
!             Update the lower triangle of A(k:n,k:n)
!
              CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                           B( K, K ), LDB, INFO )
              IF( K+KB.LE.N ) THEN
                 CALL DTRSM( 'Right', UPLO, 'Transpose', 'Non-unit', &
                             N-K-KB+1, KB, ONE, B( K, K ), LDB, &
                             A( K+KB, K ), LDA )
                 CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                             A( K, K ), LDA, B( K+KB, K ), LDB, ONE, &
                             A( K+KB, K ), LDA )
                 CALL DSYR2K( UPLO, 'No transpose', N-K-KB+1, KB, &
                              -ONE, A( K+KB, K ), LDA, B( K+KB, K ), &
                              LDB, ONE, A( K+KB, K+KB ), LDA )
                 CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                             A( K, K ), LDA, B( K+KB, K ), LDB, ONE, &
                             A( K+KB, K ), LDA )
                 CALL DTRSM( 'Left', UPLO, 'No transpose', &
                             'Non-unit', N-K-KB+1, KB, ONE, &
                             B( K+KB, K+KB ), LDB, A( K+KB, K ), &
                             LDA )
              END IF
!@CUF istat = cudaDeviceSynchronize()
   20       CONTINUE
          END IF
        ELSE
          IF( UPPER ) THEN
!
!           Compute U*A*U**T
!
            DO 30 K = 1, N, NB
              KB = MIN( N-K+1, NB )
!
!             Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
              CALL DTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', &
                          K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
              CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                          LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
              CALL DSYR2K( UPLO, 'No transpose', K-1, KB, ONE, &
                           A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, &
                           LDA )
              CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                          LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
              CALL DTRMM( 'Right', UPLO, 'Transpose', 'Non-unit', &
                          K-1, KB, ONE, B( K, K ), LDB, A( 1, K ), &
                          LDA )
!@CUF istat = cudaDeviceSynchronize()
              CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                           B( K, K ), LDB, INFO )
   30       CONTINUE
          ELSE
!
!           Compute L**T*A*L
!
            DO 40 K = 1, N, NB
              KB = MIN( N-K+1, NB )
!
!             Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
              CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', &
                          KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
              CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                          LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
              CALL DSYR2K( UPLO, 'Transpose', K-1, KB, ONE, &
                           A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, &
                           LDA )
              CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                          LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
              CALL DTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB, &
                          K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
!@CUF istat = cudaDeviceSynchronize()
              CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                           B( K, K ), LDB, INFO )
   40       CONTINUE
          END IF
        END IF
      END IF

!@CUF ISTAT = cudaMemAdvise(B, LDB*N, cudaMemAdviseUnsetReadMostly, 0)
      RETURN
!
!     End of DSYGST
!
      END
      
#ifdef _CUDA
      Integer Function cusolverConvertToEigMode(JOB) result(ijob)
      USE CUSOLVERDN
      Implicit none
      Character(1), intent(in) :: JOB
      Integer(4) :: ijob

      if (JOB.eq.'N'.or.JOB.eq.'n') then
        ijob = CUSOLVER_EIG_MODE_NOVECTOR
      else if (JOB.eq.'V'.or.JOB.eq.'v') then
        ijob = CUSOLVER_EIG_MODE_VECTOR
      end if
      end Function
!-------
      Integer Function cusolverConvertToFillMode(UPLO) result(iuplo)
      USE CUBLAS
      Implicit none
      Character(1), intent(in) :: UPLO
      Integer(4) :: iuplo

      if (UPLO.eq.'U'.or.UPLO.eq.'u') then
        iuplo = CUBLAS_FILL_MODE_UPPER
      else if (UPLO.eq.'L'.or.UPLO.eq.'l') then
        iuplo = CUBLAS_FILL_MODE_LOWER
      end if
      end Function
#endif

