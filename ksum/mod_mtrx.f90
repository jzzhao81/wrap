MODULE mod_mtrx

  INTERFACE
     SUBROUTINE mtrx_mult(arya,aryb,aryc,opaa,opbb)
       USE mod_param,ONLY : dp
       IMPLICIT NONE
       ! I/O
       COMPLEX(kind=dp), INTENT(in)  :: arya(:,:),aryb(:,:)
       COMPLEX(kind=dp), INTENT(out) :: aryc(:,:)
       ! optional I/O
       CHARACTER(len=1),INTENT(in),OPTIONAL  :: opaa, opbb
     END SUBROUTINE mtrx_mult
  END INTERFACE

END MODULE mod_mtrx


!=========>=========>=========>=========>=========>=========>=========!
! Get inverse of a complex matrix
! input dimension and matrix
! output matrix
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE mtrx_invs(ndim,mtrx)

  USE mod_param, ONLY : dp

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


! Multiplit two complex matrix by ZGEMM
SUBROUTINE mtrx_mult(arya,aryb,aryc,opaa,opbb)

  USE mod_param,ONLY : dp

  IMPLICIT NONE

  ! I/O
  COMPLEX(kind=dp), INTENT(in)  :: arya(:,:),aryb(:,:)
  COMPLEX(kind=dp), INTENT(out) :: aryc(:,:)
  ! optional I/O
  CHARACTER(len=1),INTENT(in),OPTIONAL :: opaa, opbb
  ! Work variable
  INTEGER             :: lda, ldb, ldc
  CHARACTER(len=1)    :: opa, opb
  INTEGER             :: mm, nn, kk, pp, qq
  COMPLEX(kind=dp),PARAMETER    :: alpha=CMPLX(1.d0,0.d0,dp)
  COMPLEX(kind=dp),PARAMETER    :: beta=CMPLX(0.d0,0.d0,dp)

  aryc = CMPLX(0.d0,0.d0,dp)

  ! Determin matrix operation
  IF(PRESENT(opaa)) THEN
     opa = opaa
  ELSE
     opa = 'N'
  END IF
  IF(PRESENT(opbb)) THEN
     opb = opbb
  ELSE
     opb = 'N'
  END IF

  ! Determin matrx demention
  mm = SIZE(arya, 1)
  nn = SIZE(aryb, 2)
  kk = SIZE(arya, 2)
  pp = SIZE(arya, 2)
  qq = SIZE(aryb, 1)

  lda= SIZE(arya, 1)
  ldb= SIZE(aryb, 1)
  ldc= SIZE(aryc, 1)

  IF((opa=='C').OR.(opa=='c').OR.(opa=='T').OR.(opa=='t')) THEN
     mm = SIZE(arya,2)
     kk = SIZE(arya,1)
     pp = SIZE(arya,1)
  END IF
  IF((opb=='C').OR.(opb=='c').OR.(opb=='T').OR.(opb=='t')) THEN
     nn = SIZE(aryb,1)
     qq = SIZE(aryb,2)
  END IF

  ! Check the demension parameter
  ! if 不OK, STOP
  IF((pp /= qq) .OR. (mm /= SIZE(aryc,1)) .OR. (nn /= SIZE(aryc,2)) ) THEN
     WRITE(*,*) mm, nn, kk , pp, qq
     STOP 'Error in mtrx_mult !'
  END IF

  ! call ZGEMM to multiply two arrays
  CALL zgemm(opa,opb,mm,nn,kk,alpha,arya,lda,aryb,ldb,beta,aryc,ldc)

  RETURN

END SUBROUTINE mtrx_mult
