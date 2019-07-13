MODULE mod_eig

  IMPLICIT NONE

  INTERFACE eigsys
     MODULE PROCEDURE eigsys_heig
     MODULE PROCEDURE eigsys_geig
     !module procedure eigsys_hvec
     MODULE PROCEDURE eigsys_gvec
  END INTERFACE

CONTAINS

  SUBROUTINE eigsys_heig(ndim,hams,eigs)

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.d0)

    ! I/O
    INTEGER, INTENT(in) :: ndim
    COMPLEX(kind=dp), INTENT(in) :: hams(ndim,ndim)
    REAL(kind=dp), INTENT(out) :: eigs(ndim)

    ! Work variable
    INTEGER :: lda, lwork, info
    CHARACTER :: jobz, uplo
    ! Work Array
    COMPLEX(kind=dp), ALLOCATABLE :: work(:)
    REAL(kind=dp),ALLOCATABLE ::  rwork(:)

    jobz  = 'N'; uplo  = 'U'
    lda   = ndim; lwork = 3*ndim-1

    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*ndim-2))

    eigs = 0.d0
    CALL zheev(jobz,uplo,ndim,hams,lda,eigs,work,lwork,rwork,info)

    IF(INFO.NE.0) THEN
       WRITE(*,*) 'Failure in eigsys_heig'
       WRITE(*,*) 'ZHEEV, info=',info
       STOP
    ENDIF

    DEALLOCATE(work)
    DEALLOCATE(rwork)

    RETURN

  END SUBROUTINE eigsys_heig

  SUBROUTINE eigsys_gvec(ndim,hams,eigs,evl,evr)

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.d0)
    REAL(kind=dp),PARAMETER :: dero = 0.d0
    COMPLEX(kind=dp),PARAMETER :: zero = CMPLX(dero,dero,dp)

    ! I/O
    INTEGER, INTENT(in) :: ndim
    COMPLEX(kind=dp), INTENT(in) :: hams(ndim,ndim)
    COMPLEX(kind=dp), INTENT(out) :: eigs(ndim)
    COMPLEX(kind=dp), INTENT(out) :: evl(ndim,ndim), evr(ndim,ndim)

    ! Work variable
    INTEGER   :: lda, lwork
    CHARACTER :: jobvl,jobvr
    INTEGER   :: INFO
    INTEGER   :: ldvl, ldvr
    ! Work Array
    COMPLEX(kind=dp), ALLOCATABLE :: work(:),vec(:,:)
    REAL(kind=dp),ALLOCATABLE ::  rwork(:)
    INTEGER, ALLOCATABLE :: idxarr(:)

    jobvl='V'
    jobvr='V'
    lda=ndim
    lwork=8*ndim
    ldvl=ndim
    ldvr=ndim

    ALLOCATE(work(lwork))
    ALLOCATE(rwork(8*ndim))
    ALLOCATE(vec(ndim,ndim))
    ALLOCATE(idxarr(ndim))
    work   = zero
    rwork  = dero
    vec    = hams
    idxarr = 0

    CALL zgeev(jobvl,jobvr,ndim,vec,lda,eigs,evl,ldvl,evr,ldvr,work,lwork,rwork,info)

    IF(INFO.NE.0) THEN
       WRITE(*,*) 'Failure in eigsys_gvec'
       WRITE(*,*) 'ZGEEV, info=',info
       STOP
    END IF

    ! make |vl(i)> -->  <vl(i)|
    evl = CONJG(TRANSPOSE(evl))

    ! make eigenvector orthonormal
    CALL vec_orth(evl,evr,ndim)

    ! sort eigen values according to its real part
    CALL eig_order_real_part(eigs, idxarr, ndim)
    CALL permute_eigensystem(idxarr, eigs, evl, evr, ndim)

    DEALLOCATE(work)
    DEALLOCATE(rwork)
    DEALLOCATE(vec)
    DEALLOCATE(idxarr)

    RETURN

  END SUBROUTINE eigsys_gvec


  SUBROUTINE eigsys_geig(ndim,hams,eigs)

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.d0)
    COMPLEX(kind=dp),PARAMETER :: zero=CMPLX(0.d0,0.d0,dp)

    ! I/O
    INTEGER, INTENT(in) :: ndim
    COMPLEX(kind=dp), INTENT(in) :: hams(ndim,ndim)
    COMPLEX(kind=dp), INTENT(out) :: eigs(ndim)
    INTEGER, ALLOCATABLE :: idxarr(:)
    ! Work variable
    INTEGER :: lda, lwork
    CHARACTER(len=1) :: jobvl,jobvr
    INTEGER :: INFO, ldvl, ldvr
    ! Work Array
    COMPLEX(kind=dp), ALLOCATABLE :: work(:),vec(:,:)
    REAL(kind=dp),ALLOCATABLE ::  rwork(:)
    COMPLEX(kind=dp), ALLOCATABLE :: vl(:,:), vr(:,:)

    jobvl='N'; jobvr='N'
    lda=ndim; lwork=3*ndim-1
    ldvl=1; ldvr=1

    ALLOCATE(work(lwork))
    ALLOCATE(rwork(2*ndim))
    ALLOCATE(vec(ndim,ndim))
    ALLOCATE(vl(ndim,ndim),vr(ndim,ndim))
    ALLOCATE(idxarr(ndim))

    vec = hams
    eigs = zero; work = zero
    vl = zero;   vr = zero
    rwork = 0.d0

    CALL zgeev(jobvl,jobvr,ndim,vec,lda,eigs,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

    IF(INFO.NE.0) THEN
       WRITE(*,*) 'Failure in eigsys_geig'
       WRITE(*,*) 'ZGEEV, info=',info
       STOP
    ENDIF

    ! sort eigen values according to its real part
    CALL eig_order_real_part(eigs, idxarr, ndim)

    DEALLOCATE(work)
    DEALLOCATE(vec)
    DEALLOCATE(idxarr)

    RETURN

  END SUBROUTINE eigsys_geig

  !===========================================================================
  SUBROUTINE eig_order_real_part(ev, idxarr, ndim)
    IMPLICIT NONE
!!!-----------------------------------------------------------------!!!
!!! This routine sorts complex eigenvalues of a matrix according to !!!
!!! its real parts with the smallest in the first slot and reorders !!!
!!! the matrices of left (row) and right (column) eigenvectors in a !!!
!!! corresponding manner.                                           !!!
!!!-----------------------------------------------------------------!!!
    !---------- Passed variables ----------
    COMPLEX*16, INTENT(in) :: ev(ndim)         ! Array of eigenvalues
    INTEGER, INTENT(out)   :: idxarr(ndim)     ! Index array which gives proper order
    INTEGER :: ndim                            ! Dimension of matrices
    !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
    !---------- Parameters ----------
    REAL*8, PARAMETER :: maxval = 1000.d0
    !---------- Local variables ----------
    LOGICAL, ALLOCATABLE :: sorted(:)
    REAL*8,  ALLOCATABLE :: sortonr(:)
    INTEGER :: p
    INTEGER :: q
    INTEGER :: idx
    REAL*8  :: min
    !---------- Allocate dynamic memory storage ----------
    ALLOCATE(sortonr(ndim), sorted(ndim))
    !---------- Initialize arrays ----------
    idxarr = 0
    idx    = 0
    sortonr = DBLE(ev)
    sorted = .FALSE.
    !---------- Create index array for real value ----------
    sorted = .FALSE.
    DO p = 1,ndim
       min = maxval
       DO q = 1,ndim
          IF(.NOT.sorted(q).AND.min.GT.sortonr(q)) THEN
             min = sortonr(q)
             idx = q
          ENDIF
       ENDDO
       idxarr(p) = idx
       sorted(idx) = .TRUE.
    ENDDO
    DEALLOCATE(sortonr, sorted)
    RETURN
  END SUBROUTINE eig_order_real_part

  SUBROUTINE permute_eigensystem(idxarr, ev, evl, evr, ndim)
    IMPLICIT NONE
    !---------- Passed variables ----------
    INTEGER, INTENT(in)       :: idxarr(ndim)     ! Index array which gives proper order
    COMPLEX*16, INTENT(inout) :: ev(ndim)         ! Array of eigenvalues
    COMPLEX*16, INTENT(inout) :: evl(ndim,ndim)   ! Matrix of left eigenvectors  (row)
    COMPLEX*16, INTENT(inout) :: evr(ndim,ndim)   ! Matrix of right eigenvectors (column)
    INTEGER :: ndim                               ! Dimension of matrices
    !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
    !---------- Local variables ------------------
    INTEGER :: p
    COMPLEX*16, ALLOCATABLE :: eval(:)
    COMPLEX*16, ALLOCATABLE :: evec(:,:)
    ALLOCATE(eval(ndim), evec(ndim,ndim))
    !---------- Permute the eigenvalues ----------
    DO p = 1,ndim
       eval(p) = ev(idxarr(p))
    ENDDO
    ev = eval
    !---------- Permute the right eigenvectors ----------
    DO p = 1,ndim
       evec(:,p) = evr(:,idxarr(p))
    ENDDO
    evr = evec
    !---------- Permute the left eigenvectors ----------
    DO p = 1,ndim
       evec(p,:) = evl(idxarr(p),:)
    ENDDO
    evl = evec
    !---------- Deallocate dynamic memory storage ----------
    DEALLOCATE(eval, evec)
    RETURN

  END SUBROUTINE permute_eigensystem

  ! Make vl and vr orthonormal using Gram–Schmidt process
  ! Ref: http://en.wikipedia.org/wiki/Gram–Schmidt_process
  SUBROUTINE vec_orth(vl, vr, ndim)

    IMPLICIT NONE

    INTEGER,PARAMETER :: dp = KIND(1.d0)
    COMPLEX(kind=dp),PARAMETER :: zero = CMPLX(0.d0,0.d0,dp)
    REAL(kind=dp), PARAMETER   :: eps = 1.d-5
    ! I/O
    INTEGER, INTENT(in)    ::  ndim
    COMPLEX(kind=dp), INTENT(inout) :: vr(ndim,ndim)
    COMPLEX(kind=dp), INTENT(inout) :: vl(ndim,ndim)
    ! Work variable
    INTEGER :: ii, jj
    COMPLEX(kind=dp)  :: iprod, jprod, tmp
    COMPLEX(kind=dp),EXTERNAL :: zdotu

    DO ii = 1, ndim
       DO jj = 1,ii-1
          iprod = zdotu(ndim,vl(jj,:),1,vr(:,ii),1)
          jprod = zdotu(ndim,vl(jj,:),1,vr(:,jj),1)
          IF(ABS(iprod) > eps)vr(:,ii) = vr(:,ii) - vr(:,jj)*iprod/jprod
       END DO
       DO jj = 1,ii-1
          iprod = zdotu(ndim,vl(ii,:),1,vr(:,jj),1)
          jprod = zdotu(ndim,vl(jj,:),1,vr(:,jj),1)
          IF(ABS(iprod) > eps)vl(ii,:) = vl(ii,:) - vl(jj,:)*iprod/jprod
       END DO
    END DO

    DO ii = 1, ndim
       tmp = zdotu(ndim,vl(ii,:),1,vr(:,ii),1)
       vl(ii,:) = vl(ii,:) / SQRT(tmp)
       vr(:,ii) = vr(:,ii) / SQRT(tmp)
    END DO

    RETURN
  END SUBROUTINE vec_orth


END MODULE mod_eig
