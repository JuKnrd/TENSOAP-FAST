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
   !$OMP PARALLEL
   !$OMP WORKSHARE
   prediction_lm = matmul(kte,weights) + meanval
   !$OMP END WORKSHARE
   !$OMP END PARALLEL
   prediction_lm = prediction_lm + meanval
   k = 0
   do i=1,nmol
    do mu=1,degen
     k = k + 1
     do_prediction(i,mu) = prediction_lm(k)
    enddo
   enddo

 end function

!***************************************************************************************************

 function do_scalar_kernel(PS_1,PS_2,nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0
   real*8 PS_1(nfeat_1,1,natmax_1,nmol_1),PS_2(nfeat_2,1,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8, allocatable :: PS_1_fast(:,:),PS_2_fast(:,:)
   logical hermiticity
   real*8 do_scalar_kernel(nmol_1,nmol_2)

   ! Get kernels
   if (zeta.eq.1) then
    allocate(PS_1_fast(nfeat_1,nmol_1),PS_2_fast(nfeat_2,nmol_2))
    PS_1_fast(:,:) = 0.d0
    PS_2_fast(:,:) = 0.d0
    do i=1,nmol_1
     do ii=1,int(natoms_1(i))
      PS_1_fast(:,i) = PS_1_fast(:,i) + PS_1(:,1,ii,i) / natoms_1(i)
     enddo
    enddo
    do i=1,nmol_2
     do ii=1,int(natoms_2(i))
      PS_2_fast(:,i) = PS_2_fast(:,i) + PS_2(:,1,ii,i) / natoms_2(i)
     enddo
    enddo
    do_scalar_kernel(:,:) = 0.d0
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_scalar_kernel,PS_1_fast,PS_2_fast,natoms_1,natoms_2)
     do j=j0,nmol_2
      do_scalar_kernel(i,j) = dot_product(PS_1_fast(:,i),PS_2_fast(:,j))
      if (hermiticity) do_scalar_kernel(j,i) = do_scalar_kernel(i,j)
     enddo
     !$OMP END PARALLEL DO
    enddo
    deallocate(PS_1_fast,PS_2_fast)
   else
    do_scalar_kernel(:,:) = 0.d0
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_scalar_kernel,PS_1,PS_2,natoms_1,natoms_2) PRIVATE(ii,jj)
     do j=j0,nmol_2
      do ii=1,int(natoms_1(i))
       do jj=1,int(natoms_2(j))
        do_scalar_kernel(i,j) = do_scalar_kernel(i,j) + (dot_product(PS_1(:,1,ii,i),PS_2(:,1,jj,j)))**zeta
       enddo
      enddo
      do_scalar_kernel(i,j) = do_scalar_kernel(i,j) / (natoms_1(i)*natoms_2(j))
      if (hermiticity) do_scalar_kernel(j,i) = do_scalar_kernel(i,j)
     enddo
     !$OMP END PARALLEL DO
    enddo
   endif

 end function

!***************************************************************************************************

 function do_linear_spherical_kernel(PS_1,PS_2,nmol_1,nmol_2,nfeat_1, &
     &     nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity,degen)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0,degen,mu,nu
   real*8 PS_1(nfeat_1,degen,natmax_1,nmol_1),PS_2(nfeat_2,degen,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8 do_linear_spherical_kernel(nmol_1,nmol_2,degen,degen)
   real*8, allocatable :: PS_1_fast_lm(:,:,:),PS_2_fast_lm(:,:,:)
   logical hermiticity

   do_linear_spherical_kernel(:,:,:,:) = 0.d0
   ! Fast-averaging
    allocate(PS_1_fast_lm(nfeat_1,degen,nmol_1),PS_2_fast_lm(nfeat_2,degen,nmol_2))
    PS_1_fast_lm(:,:,:) = 0.d0
    PS_2_fast_lm(:,:,:) = 0.d0
    do i=1,nmol_1
     do ii=1,int(natoms_1(i))
      PS_1_fast_lm(:,:,i) = PS_1_fast_lm(:,:,i) + PS_1(:,:,ii,i) / natoms_1(i)
     enddo
    enddo
    do i=1,nmol_2
     do ii=1,int(natoms_2(i))
      PS_2_fast_lm(:,:,i) = PS_2_fast_lm(:,:,i) + PS_2(:,:,ii,i) / natoms_2(i)
     enddo
    enddo
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_linear_spherical_kernel,PS_1_fast_lm,PS_2_fast_lm,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu)
     do j=j0,nmol_2
      do mu=1,degen
       do nu=1,degen
        do_linear_spherical_kernel(i,j,mu,nu) = (dot_product(PS_1_fast_lm(:,mu,i),PS_2_fast_lm(:,nu,j)))
        if (hermiticity) do_linear_spherical_kernel(j,i,nu,mu) = do_linear_spherical_kernel(i,j,mu,nu)
       enddo
      enddo
     enddo
     !$OMP END PARALLEL DO
    enddo

 end function

!***************************************************************************************************

 function do_nonlinear_spherical_kernel(PS_1,PS_2,PS0_1,PS0_2,nmol_1,nmol_2,nfeat_1,nfeat_2,nfeat0_1,nfeat0_2, &
     &     natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity,degen)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0,degen,mu,nu,nfeat0_1,nfeat0_2
   real*8 PS_1(nfeat_1,degen,natmax_1,nmol_1),PS_2(nfeat_2,degen,natmax_2,nmol_2)
   real*8 PS0_1(nfeat0_1,1,natmax_1,nmol_1),PS0_2(nfeat0_2,1,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2),k0
   real*8 do_nonlinear_spherical_kernel(nmol_1,nmol_2,degen,degen)
   logical hermiticity

   do_nonlinear_spherical_kernel(:,:,:,:) = 0.d0
   do i=1,nmol_1
    j0=1
    if (hermiticity) j0=i
    !$OMP PARALLEL DO SHARED(do_nonlinear_spherical_kernel,PS_1,PS_2,PS0_1,PS0_2,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu,k0)
    do j=j0,nmol_2
     do ii=1,int(natoms_1(i))
      do jj=1,int(natoms_2(j))
       k0 = (dot_product(PS0_1(:,1,ii,i),PS0_2(:,1,jj,j)))
       do mu=1,degen
        do nu=1,degen
         do_nonlinear_spherical_kernel(i,j,mu,nu) = do_nonlinear_spherical_kernel(i,j,mu,nu) + &
     &      dot_product(PS_1(:,mu,ii,i),PS_2(:,nu,jj,j)) * k0**(zeta-1)
        enddo
       enddo
      enddo
     enddo
     do_nonlinear_spherical_kernel(i,j,:,:) = do_nonlinear_spherical_kernel(i,j,:,:) / (natoms_1(i)*natoms_2(j))
     if (hermiticity) then
      do mu=1,degen
       do nu=1,degen
        do_nonlinear_spherical_kernel(j,i,nu,mu) = do_nonlinear_spherical_kernel(i,j,mu,nu)
       enddo
      enddo
     endif
    enddo
    !$OMP END PARALLEL DO
   enddo

 end function

!***************************************************************************************************
end module
