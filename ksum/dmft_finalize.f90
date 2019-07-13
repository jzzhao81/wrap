!=========>=========>=========>=========>=========>=========>=========!
! mainly deallocate global allocatable arrays
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE dmft_finalize()

  USE mod_dmft,ONLY : ek,wk,eimp,smtrx,sgm,freq,gloc,hyb,&
       &sgm_log,freq_log,hyb_log,rho,occ_dmft, sgm_old

  ! deallocate all global allatable arrays
  DEALLOCATE(ek)
  DEALLOCATE(wk)
  DEALLOCATE(eimp)
  DEALLOCATE(smtrx)
  DEALLOCATE(sgm)
  DEALLOCATE(freq)
  DEALLOCATE(gloc)
  DEALLOCATE(hyb)
  DEALLOCATE(sgm_log)
  DEALLOCATE(sgm_old)
  DEALLOCATE(freq_log)
  DEALLOCATE(hyb_log)
  DEALLOCATE(rho)
  DEALLOCATE(occ_dmft)

  RETURN

END SUBROUTINE dmft_finalize
