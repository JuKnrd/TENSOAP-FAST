module sagpr

 contains

!***************************************************************************************************

 function do_prediction(kernel,weights,meanval,degen,nmol,nenv)
  implicit none

   integer degen,nmol,nenv,i,j,k,mu,nu
   real*8 kernel(nmol,nenv,degen,degen),weights(nenv*degen),meanval,do_prediction(nmol,degen)
   real*8 kte(nmol*degen,nenv*degen),prediction_lm(nmol*degen)

   kte(:,:) = 0.d0
   do i=1,nmol
    do j=1,nenv
     do mu=1,degen
      do nu=1,degen
       kte((i-1)*degen + nu,(j-1)*degen + mu) = kernel(i,j,mu,nu)
      enddo
     enddo
    enddo
   enddo
   prediction_lm = matmul(kte,weights) + meanval
   k = 0
   do i=1,nmol
    do mu=1,degen
     k = k + 1
     do_prediction(i,mu) = prediction_lm(k)
    enddo
   enddo

 end function

!***************************************************************************************************

end module
