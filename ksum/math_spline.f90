!=====================================================================!
! interpolation subroutine with natural cubic spline method,
! warning: need more test !
!=====================================================================!

SUBROUTINE math_spline(num_in,x_in,y_in,num_out,x_out,y_out)

  USE mod_param, ONLY : dp, dzero

  IMPLICIT NONE

  INTEGER,INTENT(in) :: num_in, num_out
  REAL(kind=dp),INTENT(in) :: x_in(num_in),y_in(num_in),x_out(num_out)
  REAL(kind=dp), INTENT(out) :: y_out(num_out)
  INTEGER :: ipt

  y_out = dzero

  DO ipt=1,num_out
     CALL spline(num_in,x_in,y_in,x_out(ipt),y_out(ipt))
  END DO

  RETURN

END SUBROUTINE math_spline

SUBROUTINE spline(npt,xin,yin,xout,yout)

  USE mod_param, ONLY : dp, dzero

  IMPLICIT NONE

  INTEGER, INTENT(in) :: npt
  REAL(kind=dp),INTENT(in) :: xin(npt), yin(npt)
  REAL(kind=dp),INTENT(in) :: xout
  REAL(kind=dp),INTENT(out) :: yout

  INTEGER :: ipt
  REAL(kind=dp) :: hi(npt-1), bi(npt-1), ui(2:npt-1), vi(2:npt-1)
  REAL(kind=dp) :: y2(npt)
  REAL(kind=dp) :: a,b,c,d

  yout=dzero

  ! first we should get the second dirivation
  !===================================================================!
  DO ipt=1,npt-1
     hi(ipt)=xin(ipt+1)-xin(ipt)
     bi(ipt)=(yin(ipt+1)-yin(ipt))/hi(ipt)
  END DO

  ui(2)=2.d0*(hi(1)+hi(2))
  vi(2)=6.d0*(bi(2)-bi(1))

  DO ipt=3,npt-1
     ui(ipt)=2.d0*(hi(ipt-1)+hi(ipt))-hi(ipt-1)**2/ui(ipt-1)
     vi(ipt)=6.d0*(bi(ipt)-bi(ipt-1))-hi(ipt-1)*vi(ipt-1)/ui(ipt-1)
  END DO

  y2(1)=0.d0
  y2(npt)=0.d0

  DO ipt=npt-1,2,-1
     y2(ipt)=(vi(ipt)-hi(ipt)*y2(ipt+1))/ui(ipt)
  END DO
  !===================================================================!

  ! evaluate S
  !===================================================================!
  DO ipt=1,npt-1
     IF(xout .LE. (1.d-10+xin(ipt+1))) THEN
        EXIT
     END IF
  END DO

  a=yin(ipt)
  b=-hi(ipt)*y2(ipt+1)/6.d0-hi(ipt)*y2(ipt)/3.d0+(yin(ipt+1)-yin(ipt))&
       &/hi(ipt)
  c=y2(ipt)/2.d0
  d=(y2(ipt+1)-y2(ipt))/6.d0/hi(ipt)

  yout=a+(xout-xin(ipt))*(b+(xout-xin(ipt))*(c+(xout-xin(ipt))*d))

  RETURN

END SUBROUTINE spline
