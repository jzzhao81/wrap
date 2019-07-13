!=========>=========>=========>=========>=========>=========>=========!
! Get inverse of a complex matrix
! input dimension and matrix
! output matrix
!=========>=========>=========>=========>=========>=========>=========!

SUBROUTINE mtrx_invs(ndim,mtrx)

  USE mod_dmft, ONLY : dp
  IMPLICIT NONE

  INTEGER,INTENT(in) :: ndim
  COMPLEX(kind=dp),INTENT(inout) :: mtrx(ndim,ndim)

  INTEGER :: ipiv(ndim), info, lda, lwork
  COMPLEX(kind=dp), ALLOCATABLE :: work(:)


  lda=ndim
  lwork=ndim+1

  ALLOCATE(work(lwork))

  CALL zgetrf(ndim,ndim,mtrx,lda,ipiv,info)

  IF(info .NE. 0) WRITE(*,*) 'zgetrf info=', info

  CALL zgetri(ndim,mtrx,lda,ipiv,work,lwork,info)

  IF(info .NE. 0) WRITE(*,*) 'zgetri info=', info

  DEALLOCATE(work)

  RETURN

END SUBROUTINE mtrx_invs


SUBROUTINE mtrx_sqrt(ndim,mtrx)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ndim
  COMPLEX(kind=8), INTENT(inout) :: mtrx(ndim,ndim)
  COMPLEX(kind=8), PARAMETER :: zzero=(0.d0,0.d0)

  COMPLEX(kind=8), ALLOCATABLE :: eigvec(:,:)
  REAL(kind=8), ALLOCATABLE :: eigval(:)
  COMPLEX(kind=8), ALLOCATABLE :: mval(:,:)
  COMPLEX(kind=8), ALLOCATABLE :: rwork(:)
  COMPLEX(kind=8), ALLOCATABLE :: work(:)
  COMPLEX(kind=8), ALLOCATABLE :: tmp(:,:)
  INTEGER ::  lwork,info,ii

  lwork=2*ndim-1

  ALLOCATE(rwork(3*ndim-2))
  ALLOCATE(work(lwork))
  ALLOCATE(eigval(ndim))
  ALLOCATE(eigvec(ndim,ndim))
  ALLOCATE(tmp(ndim,ndim))
  ALLOCATE(mval(ndim,ndim))


  eigvec=mtrx

  CALL zheev('V','U',ndim,eigvec,ndim,eigval,work,lwork,rwork,info)

  mval=zzero
  DO ii=1,ndim
     mval(ii,ii)=SQRT(eigval(ii))
  END DO

  CALL zgemm('N','N',ndim,ndim,ndim,(1.d0,0.d0),eigvec,ndim,mval,ndim,(0.d0,0.d0),tmp,ndim)

  CALL zgemm('N','C',ndim,ndim,ndim,(1.d0,0.d0),tmp,ndim,eigvec,ndim,(0.d0,0.d0),mtrx,ndim)

  DEALLOCATE(rwork)
  DEALLOCATE(work)
  DEALLOCATE(eigval)
  DEALLOCATE(eigvec)
  DEALLOCATE(tmp)
  DEALLOCATE(mval)

  RETURN

END SUBROUTINE mtrx_sqrt
