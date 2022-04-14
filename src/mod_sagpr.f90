module sagpr
 use lode

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
       kte((i-1)*degen + mu,(j-1)*degen + nu) = kernel(i,j,mu,nu)
      enddo
     enddo
    enddo
   enddo

   !$OMP PARALLEL
   !$OMP WORKSHARE
   prediction_lm = matmul(kte,weights) + meanval
   !$OMP END WORKSHARE
   !$OMP END PARALLEL
   k = 0
   do i=1,nmol
    do mu=1,degen
     k = k + 1
     do_prediction(i,mu) = prediction_lm(k)
    enddo
   enddo

 end function

!***************************************************************************************************

 function do_prediction_c(kernel,weights,meanval,degen,nmol,nenv,nw)
  implicit none

   integer degen,nmol,nenv,nw,i,j,k,l,mu,nu
   real*8 kernel(nmol,nenv,degen,degen),weights(nenv*degen,nw),meanval(nw),do_prediction_c(nmol,degen,nw)
   real*8 kte(nmol*degen,nenv*degen),prediction_lm_c(nmol*degen,nw)

   kte(:,:) = 0.d0
   do i=1,nmol
    do j=1,nenv
     do mu=1,degen
      do nu=1,degen
       kte((i-1)*degen + mu,(j-1)*degen + nu) = kernel(i,j,mu,nu)
      enddo
     enddo
    enddo
   enddo

   !$OMP PARALLEL
   !$OMP WORKSHARE
   prediction_lm_c = matmul(kte,weights)
   !$OMP END WORKSHARE
   !$OMP END PARALLEL
   l = 0
   do i=1,nmol
    do mu=1,degen
     l = l + 1
     do k=1,nw
      do_prediction_c(i,mu,k) = prediction_lm_c(l,k) + meanval(k)
     enddo
    enddo
   enddo

 end function

!***************************************************************************************************
! THE FOLLOWING FUNCTION EXISTS FOR CALLING sagpr_predict
!***************************************************************************************************

 function do_prediction_rev(kernel,weights,meanval,degen,nmol,nenv)
  implicit none

   integer degen,nmol,nenv,i,j,k,mu,nu
   real*8 kernel(nmol,nenv,degen,degen),weights(nenv*degen),meanval,do_prediction_rev(nmol,degen)
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
   k = 0
   do i=1,nmol
    do mu=1,degen
     k = k + 1
     do_prediction_rev(i,mu) = prediction_lm(k)
    enddo
   enddo

 end function

!***************************************************************************************************

 function do_prediction_no_wt(kernel,meanval,degen,nmol,nenv)
  implicit none

   integer degen,nmol,nenv,i,mu
   real*8 kernel(nmol,nenv,degen,degen),meanval,do_prediction_no_wt(nmol,degen)

   do mu=1,degen
    !$OMP PARALLEL DO SHARED(kernel,meanval,do_prediction_no_wt) PRIVATE(i)
    do i=1,nmol
     do_prediction_no_wt(i,mu) = sum(kernel(i,:,mu,:)) + meanval
    enddo
    !$OMP END PARALLEL DO
   enddo

 end function

!***************************************************************************************************

 function do_scalar_kernel(PS_1,PS_2,nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0
   real*8 PS_1(nmol_1,natmax_1,1,nfeat_1),PS_2(nmol_2,natmax_2,1,nfeat_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8, allocatable :: PS_1_fast(:,:),PS_2_fast(:,:)
   logical hermiticity
   real*8 do_scalar_kernel(nmol_1,nmol_2)

   ! Get kernels
   if (zeta.eq.1) then
    allocate(PS_1_fast(nmol_1,nfeat_1),PS_2_fast(nmol_2,nfeat_2))
    PS_1_fast(:,:) = 0.d0
    PS_2_fast(:,:) = 0.d0
    do i=1,nmol_1
     do ii=1,int(natoms_1(i))
      PS_1_fast(i,:) = PS_1_fast(i,:) + PS_1(i,ii,1,:) / natoms_1(i)
     enddo
    enddo
    do i=1,nmol_2
     do ii=1,int(natoms_2(i))
      PS_2_fast(i,:) = PS_2_fast(i,:) + PS_2(i,ii,1,:) / natoms_2(i)
     enddo
    enddo
    do_scalar_kernel(:,:) = 0.d0
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_scalar_kernel,PS_1_fast,PS_2_fast,natoms_1,natoms_2)
     do j=j0,nmol_2
      do_scalar_kernel(i,j) = dot_product(PS_1_fast(i,:),PS_2_fast(j,:))
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
        do_scalar_kernel(i,j) = do_scalar_kernel(i,j) + (dot_product(PS_1(i,ii,1,:),PS_2(j,jj,1,:)))**zeta
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
   real*8 PS_1(nmol_1,natmax_1,degen,nfeat_1),PS_2(nmol_2,natmax_2,degen,nfeat_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8 do_linear_spherical_kernel(nmol_1,nmol_2,degen,degen)
   real*8, allocatable :: PS_1_fast_lm(:,:,:),PS_2_fast_lm(:,:,:)
   logical hermiticity

   do_linear_spherical_kernel(:,:,:,:) = 0.d0
   ! Fast-averaging
    allocate(PS_1_fast_lm(nmol_1,degen,nfeat_1),PS_2_fast_lm(nmol_2,degen,nfeat_2))
    PS_1_fast_lm(:,:,:) = 0.d0
    PS_2_fast_lm(:,:,:) = 0.d0
    do i=1,nmol_1
     do ii=1,int(natoms_1(i))
      PS_1_fast_lm(i,:,:) = PS_1_fast_lm(i,:,:) + PS_1(i,ii,:,:) / natoms_1(i)
     enddo
    enddo
    do i=1,nmol_2
     do ii=1,int(natoms_2(i))
      PS_2_fast_lm(i,:,:) = PS_2_fast_lm(i,:,:) + PS_2(i,ii,:,:) / natoms_2(i)
     enddo
    enddo
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_linear_spherical_kernel,PS_1_fast_lm,PS_2_fast_lm,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu)
     do j=j0,nmol_2
      do mu=1,degen
       do nu=1,degen
        do_linear_spherical_kernel(i,j,mu,nu) = (dot_product(PS_1_fast_lm(i,mu,:),PS_2_fast_lm(j,nu,:)))
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
   real*8 PS_1(nmol_1,natmax_1,degen,nfeat_1),PS_2(nmol_2,natmax_2,degen,nfeat_2)
   real*8 PS0_1(nmol_1,natmax_1,1,nfeat0_1),PS0_2(nmol_2,natmax_2,1,nfeat0_2)
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
       k0 = (dot_product(PS0_1(i,ii,1,:),PS0_2(j,jj,1,:)))
       do mu=1,degen
        do nu=1,degen
         do_nonlinear_spherical_kernel(i,j,mu,nu) = do_nonlinear_spherical_kernel(i,j,mu,nu) + &
     &      dot_product(PS_1(i,ii,mu,:),PS_2(j,jj,nu,:)) * k0**(zeta-1)
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
! THE FOLLOWING THREE FUNCTIONS ARE INCLUDED FOR POWER SPECTRA IN THE "REVERSE" ORDER, I.E.,
! THOSE THAT HAVE BEEN READ DIRECTLY FROM .BIN FILES BY SAGPR_GET_KERNEL
!***************************************************************************************************

 function do_scalar_kernel_rev(PS_1,PS_2,nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0
   real*8 PS_1(nfeat_1,1,natmax_1,nmol_1),PS_2(nfeat_2,1,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8, allocatable :: PS_1_fast(:,:),PS_2_fast(:,:)
   logical hermiticity
   real*8 do_scalar_kernel_rev(nmol_1,nmol_2)

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
    do_scalar_kernel_rev(:,:) = 0.d0
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_scalar_kernel_rev,PS_1_fast,PS_2_fast,natoms_1,natoms_2)
     do j=j0,nmol_2
      do_scalar_kernel_rev(i,j) = dot_product(PS_1_fast(:,i),PS_2_fast(:,j))
      if (hermiticity) do_scalar_kernel_rev(j,i) = do_scalar_kernel_rev(i,j)
     enddo
     !$OMP END PARALLEL DO
    enddo
    deallocate(PS_1_fast,PS_2_fast)
   else
    do_scalar_kernel_rev(:,:) = 0.d0
    do i=1,nmol_1
     j0=1
     if (hermiticity) j0=i
     !$OMP PARALLEL DO SHARED(do_scalar_kernel_rev,PS_1,PS_2,natoms_1,natoms_2) PRIVATE(ii,jj)
     do j=j0,nmol_2
      do ii=1,int(natoms_1(i))
       do jj=1,int(natoms_2(j))
        do_scalar_kernel_rev(i,j) = do_scalar_kernel_rev(i,j) + (dot_product(PS_1(:,1,ii,i),PS_2(:,1,jj,j)))**zeta
       enddo
      enddo
      do_scalar_kernel_rev(i,j) = do_scalar_kernel_rev(i,j) / (natoms_1(i)*natoms_2(j))
      if (hermiticity) do_scalar_kernel_rev(j,i) = do_scalar_kernel_rev(i,j)
     enddo
     !$OMP END PARALLEL DO
    enddo
   endif

 end function

!***************************************************************************************************

 function do_linear_spherical_kernel_rev(PS_1,PS_2,nmol_1,nmol_2,nfeat_1, &
     &     nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity,degen)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0,degen,mu,nu
   real*8 PS_1(nfeat_1,degen,natmax_1,nmol_1),PS_2(nfeat_2,degen,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2)
   real*8 do_linear_spherical_kernel_rev(nmol_1,nmol_2,degen,degen)
   real*8, allocatable :: PS_1_fast_lm(:,:,:),PS_2_fast_lm(:,:,:)
   logical hermiticity

   do_linear_spherical_kernel_rev(:,:,:,:) = 0.d0
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
     !$OMP PARALLEL DO SHARED(do_linear_spherical_kernel_rev,PS_1_fast_lm,PS_2_fast_lm,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu)
     do j=j0,nmol_2
      do mu=1,degen
       do nu=1,degen
        do_linear_spherical_kernel_rev(i,j,mu,nu) = (dot_product(PS_1_fast_lm(:,mu,i),PS_2_fast_lm(:,nu,j)))
        if (hermiticity) do_linear_spherical_kernel_rev(j,i,nu,mu) = do_linear_spherical_kernel_rev(i,j,mu,nu)
       enddo
      enddo
     enddo
     !$OMP END PARALLEL DO
    enddo

 end function

!***************************************************************************************************

 function do_nonlinear_spherical_kernel_rev(PS_1,PS_2,PS0_1,PS0_2,nmol_1,nmol_2,nfeat_1,nfeat_2,nfeat0_1,nfeat0_2, &
     &     natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity,degen)
  implicit none

   integer nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,zeta,i,j,ii,jj,j0,degen,mu,nu,nfeat0_1,nfeat0_2
   real*8 PS_1(nfeat_1,degen,natmax_1,nmol_1),PS_2(nfeat_2,degen,natmax_2,nmol_2)
   real*8 PS0_1(nfeat0_1,1,natmax_1,nmol_1),PS0_2(nfeat0_2,1,natmax_2,nmol_2)
   real*8 natoms_1(nmol_1),natoms_2(nmol_2),k0
   real*8 do_nonlinear_spherical_kernel_rev(nmol_1,nmol_2,degen,degen)
   logical hermiticity

   do_nonlinear_spherical_kernel_rev(:,:,:,:) = 0.d0
   do i=1,nmol_1
    j0=1
    if (hermiticity) j0=i
    !$OMP PARALLEL DO SHARED(do_nonlinear_spherical_kernel_rev,PS_1,PS_2,PS0_1,PS0_2,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu,k0)
    do j=j0,nmol_2
     do ii=1,int(natoms_1(i))
      do jj=1,int(natoms_2(j))
       k0 = (dot_product(PS0_1(:,1,ii,i),PS0_2(:,1,jj,j)))
       do mu=1,degen
        do nu=1,degen
         do_nonlinear_spherical_kernel_rev(i,j,mu,nu) = do_nonlinear_spherical_kernel_rev(i,j,mu,nu) + &
     &      dot_product(PS_1(:,mu,ii,i),PS_2(:,nu,jj,j)) * k0**(zeta-1)
        enddo
       enddo
      enddo
     enddo
     do_nonlinear_spherical_kernel_rev(i,j,:,:) = do_nonlinear_spherical_kernel_rev(i,j,:,:) / (natoms_1(i)*natoms_2(j))
     if (hermiticity) then
      do mu=1,degen
       do nu=1,degen
        do_nonlinear_spherical_kernel_rev(j,i,nu,mu) = do_nonlinear_spherical_kernel_rev(i,j,mu,nu)
       enddo
      enddo
     endif
    enddo
    !$OMP END PARALLEL DO
   enddo

 end function

!***************************************************************************************************

 subroutine do_power_spectrum(PS,components,xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,mult_by_A,w3j,all_species,all_centres,isLODE,LODE_params)
  use SHTOOLS, only: wigner3j
  implicit none

  complex*16, allocatable, intent(inout) :: PS(:,:,:,:)
  integer, allocatable, intent(inout) :: components(:,:)
  real*8, allocatable, intent(inout) :: w3j(:,:,:,:)
  integer nframes,natmax,lm,nmax,lmax,ncut,degen,featsize,nnmax,nsmax,ispe,i,j,k,nspecies,n1,n2,l1,l2,l,mu,nu
  integer llmax,ps_shape(4),m,n,nn,jmin,jmax,im,ia,ib,mm
  real*8 xyz(nframes,natmax,3),rs(3),rcut,sg,cell(nframes,3,3)
  complex*16 sparsification(2,ncut,ncut)
  real*8 sigma(nmax),overlap(nmax,nmax),eigenval(nmax),diagsqrt(nmax,nmax),orthomatrix(nmax,nmax),inner
  character(len=4) atname(nframes,natmax)
  integer natoms(nframes),ipiv(nmax)
  logical periodic,mult_by_A
  integer, allocatable :: all_indices(:,:,:),nneighmax(:,:),ncen(:),lvalues(:,:)
  integer, parameter :: lwmax = 10000
  integer info,lwork,work(lwmax)
  complex*16, allocatable :: omega1(:,:,:,:,:),omega2(:,:,:,:,:),harmonic(:,:,:,:,:),omegatrue(:,:,:,:,:)
  complex*16, allocatable :: omegatrue2(:,:,:,:,:),ps_row(:,:,:)
  complex*16, allocatable :: ot1(:,:),ot2(:,:)
  complex*16, allocatable :: CC(:,:),inner_mu(:,:)
  real*8, allocatable :: orthoradint(:,:,:,:,:),tmp_3j(:)
  integer, allocatable :: index_list(:)
  complex*16, allocatable :: om(:,:),ch(:,:)
  real*8, allocatable :: ww(:,:)
  integer cr,ts,tf
  real*8 rate,dr(3),invcell(3,3),sv(3)
  logical all_species(nelements),all_centres(nelements),isLODE
  type(LODE_Model), intent(in) :: LODE_params

  call system_clock(count_rate=cr)
  rate = real(cr)

  ! Get maximum number of neighbours
  if (.not.periodic) then
   nnmax = natmax
  else
   nnmax = 0
   do i=1,nframes
    ! Invert unit cell
    invcell(:,:) = cell(i,:,:)
    lwork = 1000
    call DGETRF(3,3,invcell,3,ipiv,info)
    call DGETRI(3,invcell,3,ipiv,work,lwork,info)
    do j=1,natmax
     nn = 0
     do k=1,natmax
      if (k.ne.j) then
       dr(:) = xyz(i,j,:) - xyz(i,k,:)
       sv = matmul(invcell,dr)
       do l=1,3
        sv(l) = sv(l) - nint(sv(l))
       enddo
       dr = matmul(cell(i,:,:),sv)
       if (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) .le. rcut*rcut) nn = nn + 1
      endif
     enddo
     nnmax = max(nnmax,nn)
    enddo
   enddo
  endif

  ! List indices for atoms of the same species
  nsmax = nelements
  allocate(all_indices(nframes,nsmax,natmax))
  all_indices(:,:,:) = 0
  do i=1,nframes
   do j=1,natmax
    do ispe=1,nelements
     if (trim(adjustl(atname(i,j))).eq.trim(adjustl(atomic_names(ispe)))) then
      do k=1,natmax
       if (all_indices(i,ispe,k).eq.0) then 
        all_indices(i,ispe,k) = j
        exit
       endif
      enddo
     endif
    enddo
   enddo
  enddo

  ! Maximum number of nearest neighbours
  allocate(nneighmax(nframes,nsmax))
  nneighmax(:,:) = 0
  do i=1,nframes
   do j=1,nsmax
    do k=1,natmax
     if (all_indices(i,j,k).gt.0) nneighmax(i,j) = nneighmax(i,j) + 1
    enddo
   enddo
  enddo

  ! Get list of centres and species
  if (count(all_centres).eq.0) then
   do i=1,nframes
    do j=1,natmax
     do k=1,nelements
      if (trim(adjustl(atname(i,j))).eq.atomic_names(k)) all_centres(k) = .true.
     enddo
    enddo
   enddo
  endif
  allocate(ncen(nframes))
  ncen(:) = 0
  do i=1,nframes
   do j=1,nelements
    if (all_centres(j)) then
!     do k=1,natmax
!      !if (atn
!     enddo
     ncen(i) = ncen(i) + count(atname(i,:).eq.atomic_names(j))
    endif
   enddo
  enddo

  if (count(all_species).eq.0) then
   do i=1,nframes
    do j=1,natmax
     do k=1,nelements
      if (trim(adjustl(atname(i,j))).eq.atomic_names(k)) all_species(k) = .true.
     enddo
    enddo
   enddo
  endif
  nspecies = count(all_species)

  ! Set up orthogonal matrix
  sigma(:) = 0.d0
  do i=1,nmax
   sigma(i) = max(sqrt(dble(i-1)),1.0) * rcut / dble(nmax)
  enddo
  do i=1,nmax
   do j=1,nmax
    n1=i-1
    n2=j-1
    overlap(i,j) = (0.5d0 / (sigma(i))**2 + 0.5/(sigma(j))**2)**(-0.5d0*(3.d0 + n1 + n2)) &
     &     /(sigma(i)**n1 * sigma(j)**n2)* gamma(0.5d0*(3.d0 + n1 + n2))/ &
     &     ((sigma(i)*sigma(j))**1.5d0 * sqrt(gamma(1.5d0 + n1)*gamma(1.5d0 + n2)))
   enddo
  enddo

  lwork = lwmax
  ! Get eigenvalues and eigenvectors
  call DSYEV('Vectors','Upper',nmax,overlap,nmax,eigenval,work,lwork,info)
  diagsqrt(:,:) = 0.d0
  do i=1,nmax
   diagsqrt(i,i) = sqrt(eigenval(i))
  enddo
  orthomatrix = matmul(overlap,matmul(diagsqrt,transpose(overlap)))
  ! Invert square root of overlap matrix to get orthogonal matrix
  call DGETRF(nmax,nmax,orthomatrix,nmax,ipiv,info)
  call DGETRI(nmax,orthomatrix,nmax,ipiv,work,lwork,info)

  ! Set up arrays for time-saving in spherical power spectra
  if (lm.gt.0) then
   llmax = 0
   do l1=0,lmax
    do l2=0,lmax
     if (mod((lm + l1 + l2),2).eq.0) then
      if (abs(l2-lm).le.l1 .and. l1.le.(l2+lm)) llmax = llmax + 1
     endif
    enddo
   enddo
   allocate(lvalues(llmax,2))
   llmax = 0
   do l1=0,lmax
    do l2=0,lmax
     if (mod((lm + l1 + l2),2).eq.0) then
      if (abs(l2-lm).le.l1 .and. l1.le.(l2+lm)) then
       llmax = llmax + 1
       lvalues(llmax,:) = (/l1,l2/)
      endif
     endif
    enddo
   enddo
  endif

  ! Allocate power spectrum array
  degen = 2*lm + 1
  if (ncut.ne.-1) then
   featsize = ncut
  else
   if (lm.eq.0) then
    featsize = nspecies*nspecies*nmax**2*(lmax+1)
   else
    featsize = nspecies*nspecies*nmax**2*llmax
   endif
  endif
  if (allocated(PS)) deallocate(PS)
  allocate(PS(nframes,natmax,degen,featsize))
  PS(:,:,:,:) = (0.d0,0.d0)

  ! Get list of components
  if (.not. allocated(components) .and. ncut.gt.0) then
   if (lm.eq.0) then
    allocate(components(ncut,5))
    n = 0
    do i=1,nspecies
     do j=1,nspecies
      do k=1,nmax
       do l=1,nmax
        do m=1,lmax+1
         n = n + 1
         ! Check if this component is required
         do nn=1,ncut
          ! Remember we have to add 1 because we're dealing with python arrays
          if (sparsification(1,nn,1)+1.eq.n) then
           components(nn,:) = (/i,j,k,l,m/)
          endif
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   else
    allocate(components(ncut,6))
    n = 0
    do i=1,nspecies
     do j=1,nspecies
      do k=1,nmax
       do l=1,nmax
        do m=1,llmax
         n = n + 1
         ! Check if this component is required
         do nn=1,ncut
          ! Remember we have to add 1 because we're dealing with python arrays
          if (sparsification(1,nn,1)+1.eq.n) then
           components(nn,:) = (/i,j,k,l,lvalues(m,1),lvalues(m,2)/)
          endif
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   endif
  endif

  ! If necessary, pre-compute the Wigner 3j symbols
  if (lm.gt.0) then
   if (.not. allocated(w3j)) then
    allocate(w3j(2*lm+1,lmax+1,lmax+1,2*lmax+1))
    w3j(:,:,:,:) = 0.d0
    do l1=0,lmax
     do l2=0,lmax
      do m=0,2*l1
       do mu=0,2*lm
        if (allocated(tmp_3j)) deallocate(tmp_3j)
        allocate(tmp_3j(l2+l1+1))
        tmp_3j(:) = 0.d0
        call wigner3j(tmp_3j,jmin,jmax,l2,l1,mu-lm,m-l1-mu+lm,-m+l1)
        if (lm.ge.jmin .and. lm.le.jmax) w3j(mu+1,l1+1,l2+1,m+1) = tmp_3j(lm-jmin+1) * (-1.d0)**(m-l1)
       enddo
      enddo
     enddo
    enddo
   endif
  endif

  if (isLODE .and. periodic) then
   ! If necessary, set up or recalculate the G-vectors
   call get_Gvectors(LODE_params,cell,invcell,nmax,lmax,rcut,sigma,orthomatrix)
  endif

  ! Do the power spectrum computation
  do i=1,nframes

   ! Get omega, harmonic and radint matrices
   allocate(omega1(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),harmonic(natoms(i),nspecies,lmax+1,2*lmax+1,nnmax),&
     &     orthoradint(natoms(i),nspecies,lmax+1,nmax,nnmax),omega2(natoms(i),nspecies,nmax,lmax+1,2*lmax+1))
   call initsoap(omega1,harmonic,orthoradint,natoms(i),nspecies,nmax,lmax,nnmax,periodic,all_indices(i,:,:), &
     &     nneighmax(i,:),natmax,nsmax,cell(i,:,:),rs,sg,rcut,xyz(i,:,:),sigma,orthomatrix,all_species,all_centres)

   if (.not.isLODE) then
    omega2 = omega1
   else
    if (.not.periodic) then
     call direct_potential(omega2,natoms(i),nspecies,nmax,lmax,nnmax,all_indices(i,:,:), &
     &     nneighmax(i,:),natmax,nsmax,sg,rcut,xyz(i,:,:),sigma,orthomatrix,all_species, &
     &     all_centres,LODE_params%radsize,LODE_params%lebsize)
    else
     stop 'PERIODIC LODE NOT YET IMPLEMENTED!'
    endif
   endif

   ! Compute power spectrum
   if (lm.eq.0) then
    ! Scalar
    allocate(omegatrue(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),omegatrue2(natoms(i),nspecies,nmax,lmax+1,2*lmax+1))
    do l=0,lmax
     omegatrue(:,:,:,l+1,:)  = omega2(:,:,:,l+1,:) / dsqrt(dsqrt(2.d0*l+1.d0))
     omegatrue2(:,:,:,l+1,:) = omega1(:,:,:,l+1,:) / dsqrt(dsqrt(2.d0*l+1.d0))
    enddo
    if (ncut.gt.0) then
     !$OMP PARALLEL DO SHARED(components,PS,omegatrue,omegatrue2,ncut,natoms,i) PRIVATE(j,k,ot1,ot2)
     do j=1,ncut
      allocate(ot1(natoms(i),2*lmax+1),ot2(natoms(i),2*lmax+1))
      ot1(:,:) = omegatrue(:,components(j,1),components(j,3),components(j,5),:)
      ot2(:,:) = omegatrue2(:,components(j,2),components(j,4),components(j,5),:)
      do k=1,natoms(i)
       PS(i,k,1,j) = dot_product(ot1(k,:),ot2(k,:))
      enddo
      deallocate(ot1,ot2)
     enddo
     !$OMP END PARALLEL DO
    else
     stop 'ERROR: no sparsification information given; this is not recommended!'
    endif
    deallocate(omegatrue,omegatrue2)

   else
     ! Spherical

     if (ncut.gt.0) then
      !$OMP PARALLEL DO SHARED(components,PS,omega1,omega2,orthoradint,harmonic,w3j) PRIVATE(j,ia,ib,nn,mm,l1,l2,k,l,im,n,om,ww,ch)
      do j=1,ncut
       allocate(om(natoms(i),2*lmax+1),ww(2*lm+1,2*lmax+1),ch(natoms(i),2*lmax+1))
       ia = components(j,1)
       ib = components(j,2)
       nn = components(j,3)
       mm = components(j,4)
       l1 = components(j,5)
       l2 = components(j,6)
       ww(:,:)  = w3j(:,l1+1,l2+1,:)
       om(:,:)  = omega2(:,ia,nn,l1+1,:)
       do l=1,natoms(i)
        do im=1,2*lmax+1
         ch(l,im) = dot_product(harmonic(l,ib,l2+1,im,:),orthoradint(l,ib,l2+1,mm,:))
        enddo
       enddo
       do k=1,2*lm+1
        do l=1,natoms(i)
         do im=1,2*lmax+1
          if (abs(im-l1-k+lm).le.l2) then
            PS(i,l,k,j) = PS(i,l,k,j) + (om(l,im)*ch(l,l2+im-l1-k+lm+1)*ww(k,im))
          endif
         enddo
        enddo
       enddo
       deallocate(om,ww,ch)
      enddo
      !$OMP END PARALLEL DO
     else
      stop 'ERROR: no sparsification information given; this is not recommended!'
     endif
   endif

   ! Deallocate
   deallocate(omega1,harmonic,orthoradint,omega2)

  enddo

  ! Multiply by A matrix
  if (mult_by_A) then
   do i=1,nframes
    !$OMP PARALLEL DO SHARED(PS,sparsification) PRIVATE(j,k)
    do j=1,natmax
     do k=1,degen
      PS(i,j,k,:) = matmul(PS(i,j,k,:),sparsification(2,:,:))
     enddo
    enddo
    !$OMP END PARALLEL DO
   enddo
  endif

  ! Make power spectrum real and normalize it
  if (lm.eq.0) then
   ! Scalar
   PS(:,:,:,:) = real(PS(:,:,:,:))
   do i=1,nframes
    do j=1,natoms(i)
     inner = dot_product(PS(i,j,1,:),PS(i,j,1,:))
     PS(i,j,1,:) = PS(i,j,1,:) / dsqrt(inner)
    enddo
   enddo
  else
   ! Spherical
   ! Get complex to real transformation
   allocate(CC(degen,degen),inner_mu(degen,degen))
   CC = complex_to_real_matrix(degen)
   ! Apply transformation
   do i=1,nframes
    do j=1,natoms(i)
     do k=1,ncut
      PS(i,j,:,k) = real(matmul(CC,PS(i,j,:,k)))
     enddo
    enddo
   enddo
   do i=1,nframes
    do j=1,natoms(i)
     inner = 0.d0
     do mu=1,degen
      do nu=1,degen
       inner_mu(mu,nu) = dot_product(PS(i,j,mu,:),PS(i,j,nu,:))
       inner = inner + abs(inner_mu(mu,nu))
      enddo
     enddo
     PS(i,j,:,:) = PS(i,j,:,:) / dsqrt(norm2(real(inner_mu)))
    enddo
   enddo
   deallocate(CC,inner_mu)
  endif

  ! Reorder the power spectrum so that the ordering of atoms matches their positions in the frame
  allocate(ps_row(natmax,degen,featsize),index_list(natmax))
  do i=1,nframes
   ps_row(:,:,:) = 0.d0
   index_list(:) = 0
   l = 0
   do j=1,nelements
    do k=1,natmax
     if (all_indices(i,j,k).gt.0) then
      l = l + 1
      index_list(l) = all_indices(i,j,k)
     endif
    enddo
   enddo
   do j=1,natoms(i)
    ps_row(index_list(j),:,:) = PS(i,j,:,:)
   enddo
   PS(i,:,:,:) = ps_row(:,:,:)
  enddo
  deallocate(ps_row,index_list)

  deallocate(all_indices,nneighmax,ncen)
  if (allocated(lvalues)) deallocate(lvalues)
  if (allocated(tmp_3j)) deallocate(tmp_3j)

 end subroutine

!***************************************************************************************************
! The subroutine below is for building the power spectrum that goes with zeta>1
!***************************************************************************************************

 subroutine do_power_spectrum_scalar(PS0,components0,xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,mult_by_A,all_species,all_centres,isLODE,LODE_params)
  use SHTOOLS, only: wigner3j
  implicit none

  complex*16, allocatable, intent(inout) :: PS0(:,:,:,:)
  integer, allocatable, intent(inout) :: components0(:,:)
  integer nframes,natmax,lm,nmax,lmax,ncut,degen,featsize,nnmax,nsmax,ispe,i,j,k,nspecies,n1,n2,l
  integer llmax,ps_shape(4),m,n,nn
  real*8 xyz(nframes,natmax,3),rs(3),rcut,sg,cell(nframes,3,3)
  complex*16 sparsification(2,ncut,ncut)
  real*8 sigma(nmax),overlap(nmax,nmax),eigenval(nmax),diagsqrt(nmax,nmax),orthomatrix(nmax,nmax),inner
  character(len=4) atname(nframes,natmax)
  integer natoms(nframes),ipiv(nmax)
  logical periodic,mult_by_A
  integer, allocatable :: all_indices(:,:,:),nneighmax(:,:),ncen(:)
  integer, parameter :: lwmax = 10000
  integer info,lwork,work(lwmax)
  complex*16, allocatable :: omega1(:,:,:,:,:),omega2(:,:,:,:,:),harmonic(:,:,:,:,:),omegatrue(:,:,:,:,:)
  complex*16, allocatable :: omegatrue2(:,:,:,:,:),ps_row(:,:,:)
  complex*16, allocatable :: ot1(:,:),ot2(:,:)
  real*8, allocatable :: orthoradint(:,:,:,:,:),tmp_3j(:)
  integer, allocatable :: index_list(:)
  real*8 dr(3),sv(3),invcell(3,3)
  logical all_species(nelements),all_centres(nelements),isLODE
  type(LODE_Model), intent(in) :: LODE_params

  if (lm.ne.0) stop 'ERROR: scalar power spectrum has been called with non-scalar argument!'

  ! Get maximum number of neighbours
  if (.not.periodic) then
   nnmax = natmax
  else
   nnmax = 0
   do i=1,nframes
    ! Invert unit cell
    invcell(:,:) = cell(i,:,:)
    lwork = 1000
    call DGETRF(3,3,invcell,3,ipiv,info)
    call DGETRI(3,invcell,3,ipiv,work,lwork,info)
    do j=1,natmax
     nn = 0
     do k=1,natmax
      if (k.ne.j) then
       dr(:) = xyz(i,j,:) - xyz(i,k,:)
       sv = matmul(invcell,dr)
       do l=1,3
        sv(l) = sv(l) - nint(sv(l))
       enddo
       dr = matmul(cell(i,:,:),sv)
       if (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) .le. rcut*rcut) nn = nn + 1
      endif
     enddo
     nnmax = max(nnmax,nn)
    enddo
   enddo
  endif

  ! List indices for atoms of the same species
  nsmax = nelements
  allocate(all_indices(nframes,nsmax,natmax))
  all_indices(:,:,:) = 0
  do i=1,nframes
   do j=1,natmax
    do ispe=1,nelements
     if (trim(adjustl(atname(i,j))).eq.trim(adjustl(atomic_names(ispe)))) then
      do k=1,natmax
       if (all_indices(i,ispe,k).eq.0) then 
        all_indices(i,ispe,k) = j
        exit
       endif
      enddo
     endif
    enddo
   enddo
  enddo

  ! Maximum number of nearest neighbours
  allocate(nneighmax(nframes,nsmax))
  nneighmax(:,:) = 0
  do i=1,nframes
   do j=1,nsmax
    do k=1,natmax
     if (all_indices(i,j,k).gt.0) nneighmax(i,j) = nneighmax(i,j) + 1
    enddo
   enddo
  enddo

  ! Get list of centres and species
  if (count(all_centres).eq.0) then
   do i=1,nframes
    do j=1,natmax
     do k=1,nelements
      if (trim(adjustl(atname(i,j))).eq.atomic_names(k)) all_centres(k) = .true.
     enddo
    enddo
   enddo
  endif
  allocate(ncen(nframes))
  ncen(:) = 0
  do i=1,nframes
   do j=1,nelements
    if (all_centres(j)) then
!     do k=1,natmax
!      !if (atn
!     enddo
     ncen(i) = ncen(i) + count(atname(i,:).eq.atomic_names(j))
    endif
   enddo
  enddo

  if (count(all_species).eq.0) then
   do i=1,nframes
    do j=1,natmax
     do k=1,nelements
      if (trim(adjustl(atname(i,j))).eq.atomic_names(k)) all_species(k) = .true.
     enddo
    enddo
   enddo
  endif
  nspecies = count(all_species)

  ! Set up orthogonal matrix
  sigma(:) = 0.d0
  do i=1,nmax
   sigma(i) = max(sqrt(dble(i-1)),1.0) * rcut / dble(nmax)
  enddo
  do i=1,nmax
   do j=1,nmax
    n1=i-1
    n2=j-1
    overlap(i,j) = (0.5d0 / (sigma(i))**2 + 0.5/(sigma(j))**2)**(-0.5d0*(3.d0 + n1 + n2)) &
     &     /(sigma(i)**n1 * sigma(j)**n2)* gamma(0.5d0*(3.d0 + n1 + n2))/ &
     &     ((sigma(i)*sigma(j))**1.5d0 * sqrt(gamma(1.5d0 + n1)*gamma(1.5d0 + n2)))
   enddo
  enddo

  lwork = lwmax
  ! Get eigenvalues and eigenvectors
  call DSYEV('Vectors','Upper',nmax,overlap,nmax,eigenval,work,lwork,info)
  diagsqrt(:,:) = 0.d0
  do i=1,nmax
   diagsqrt(i,i) = sqrt(eigenval(i))
  enddo
  orthomatrix = matmul(overlap,matmul(diagsqrt,transpose(overlap)))
  ! Invert square root of overlap matrix to get orthogonal matrix
  call DGETRF(nmax,nmax,orthomatrix,nmax,ipiv,info)
  call DGETRI(nmax,orthomatrix,nmax,ipiv,work,lwork,info)

  ! Allocate power spectrum array
  degen = 2*lm + 1
  if (ncut.ne.-1) then
   featsize = ncut
  else
   if (lm.eq.0) then
    featsize = nspecies*nspecies*nmax**2*(lmax+1)
   else
    featsize = nspecies*nspecies*nmax**2*llmax
   endif
  endif
  if (allocated(PS0)) deallocate(PS0)
  allocate(PS0(nframes,natmax,degen,featsize))
  PS0(:,:,:,:) = (0.d0,0.d0)

  ! Get list of components
  if (.not. allocated(components0) .and. ncut.gt.0) then
   allocate(components0(ncut,5))
   n = 0
   do i=1,nspecies
    do j=1,nspecies
    do k=1,nmax
      do l=1,nmax
       do m=1,lmax+1
        n = n + 1
        ! Check if this component is required
        do nn=1,ncut
         ! Remember we have to add 1 because we're dealing with python arrays
         if (sparsification(1,nn,1)+1.eq.n) then
          components0(nn,:) = (/i,j,k,l,m/)
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  endif

  if (isLODE .and. periodic) then
   ! If necessary, set up or recalculate the G-vectors
   call get_Gvectors(LODE_params,cell,invcell,nmax,lmax,rcut,sigma,orthomatrix)
  endif

  ! Do the power spectrum computation
  do i=1,nframes

   ! Get omega, harmonic and radint matrices
   allocate(omega1(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),harmonic(natoms(i),nspecies,lmax+1,2*lmax+1,nnmax),&
     &     orthoradint(natoms(i),nspecies,lmax+1,nmax,nnmax),omega2(natoms(i),nspecies,nmax,lmax+1,2*lmax+1))
   call initsoap(omega1,harmonic,orthoradint,natoms(i),nspecies,nmax,lmax,nnmax,periodic,all_indices(i,:,:), &
     &     nneighmax(i,:),natmax,nsmax,cell(i,:,:),rs,sg,rcut,xyz(i,:,:),sigma,orthomatrix,all_species,all_centres)

   if (.not.isLODE) then
    omega2 = omega1
   else
    if (.not.periodic) then
     call direct_potential(omega2,natoms(i),nspecies,nmax,lmax,nnmax,all_indices(i,:,:), &
     &     nneighmax(i,:),natmax,nsmax,sg,rcut,xyz(i,:,:),sigma,orthomatrix,all_species, &
     &     all_centres,LODE_params%radsize,LODE_params%lebsize)
    else
     stop 'PERIODIC LODE NOT YET IMPLEMENTED!'
    endif
   endif

   ! Compute power spectrum
   allocate(omegatrue(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),omegatrue2(natoms(i),nspecies,nmax,lmax+1,2*lmax+1))
   do l=0,lmax
    omegatrue(:,:,:,l+1,:)  = omega2(:,:,:,l+1,:) / dsqrt(dsqrt(2.d0*l+1.d0))
    omegatrue2(:,:,:,l+1,:) = omega1(:,:,:,l+1,:) / dsqrt(dsqrt(2.d0*l+1.d0))
   enddo
   if (ncut.gt.0) then
    !$OMP PARALLEL DO SHARED(components0,PS0,omegatrue,omegatrue2,ncut,natoms,i) PRIVATE(j,k,ot1,ot2)
    do j=1,ncut
     allocate(ot1(natoms(i),2*lmax+1),ot2(natoms(i),2*lmax+1))
     ot1(:,:) = omegatrue(:,components0(j,1),components0(j,3),components0(j,5),:)
     ot2(:,:) = omegatrue2(:,components0(j,2),components0(j,4),components0(j,5),:)
     do k=1,natoms(i)
      PS0(i,k,1,j) = dot_product(ot1(k,:),ot2(k,:))
     enddo
     deallocate(ot1,ot2)
    enddo
    !$OMP END PARALLEL DO
   else
    stop 'ERROR: no sparsification information given; this is not recommended!'
   endif
   deallocate(omegatrue,omegatrue2)

   ! Deallocate
   deallocate(omega1,harmonic,orthoradint,omega2)

  enddo

  ! Multiply by A matrix
  if (mult_by_A) then
   do i=1,nframes
    !$OMP PARALLEL DO SHARED(PS0,sparsification) PRIVATE(j,k)
    do j=1,natmax
     do k=1,degen
      PS0(i,j,k,:) = matmul(PS0(i,j,k,:),sparsification(2,:,:))
     enddo
    enddo
    !$OMP END PARALLEL DO
   enddo
  endif

  ! Make power spectrum real and normalize it
  PS0(:,:,:,:) = real(PS0(:,:,:,:))
  do i=1,nframes
   do j=1,natoms(i)
    inner = dot_product(PS0(i,j,1,:),PS0(i,j,1,:))
    PS0(i,j,1,:) = PS0(i,j,1,:) / dsqrt(inner)
   enddo
  enddo

  ! Reorder the power spectrum so that the ordering of atoms matches their positions in the frame
  allocate(ps_row(natmax,degen,featsize),index_list(natmax))
  do i=1,nframes
   ps_row(:,:,:) = 0.d0
   index_list(:) = 0
   l = 0
   do j=1,nelements
    do k=1,natmax
     if (all_indices(i,j,k).gt.0) then
      l = l + 1
      index_list(l) = all_indices(i,j,k)
     endif
    enddo
   enddo
   do j=1,natoms(i)
    ps_row(index_list(j),:,:) = PS0(i,j,:,:)
   enddo
   PS0(i,:,:,:) = ps_row(:,:,:)
  enddo
  deallocate(ps_row,index_list)

  deallocate(all_indices,nneighmax,ncen)
  if (allocated(tmp_3j)) deallocate(tmp_3j)

 end subroutine

!***************************************************************************************************

 subroutine initsoap(omega,harmonic,orthoradint,natoms,nspecies,nmax,lmax,nnmax,periodic,all_indices,nneighmax, &
     &     natmax,nsmax,cell,rs,sg,rcut,xyz,sigma,orthomatrix,all_species,all_centres)
  implicit none

   integer natoms,nspecies,nmax,lmax,nnmax,natmax,nsmax,iat,ncentype,icentype,icen,cen,n,ispe
   integer ineigh,neigh,lval,im,mval,n1,n2,l,i,k,nn,ncell,ia,ib,ic
   real*8 rcut2,rx,ry,rz,r2,rdist,cth,ph,normfact,sigmafact,rv(3),sv(3),rcv(3)
   complex*16 omega(natoms,nspecies,nmax,lmax+1,2*lmax+1),harmonic(natoms,nspecies,lmax+1,2*lmax+1,nnmax)
   real*8 orthoradint(natoms,nspecies,lmax+1,nmax,nnmax),f1
   real*8 cell(3,3),rs(3),sg,rcut,xyz(natmax,3)
   real*8, allocatable :: radint(:,:,:,:,:),efact(:,:,:),length(:,:,:)
   real*8 sigma(nmax),orthomatrix(nmax,nmax),alpha,sg2,radial_c,radial_r0,radial_m,invcell(3,3)
   integer, allocatable :: nneigh(:,:)
   integer all_indices(nsmax,natmax),nneighmax(nsmax),ipiv(3),info,lwork,work(1000)
   logical periodic
   integer ts,tf,cr
   real*8 rate
   logical all_species(nelements),all_centres(nelements)

   sg2 = sg*sg
   alpha = 1.d0 / (2.d0 * sg2)
   radial_c = rs(1)
   radial_r0 = rs(2)
   radial_m = rs(3)
   rcut2 = rcut*rcut
   ncell = 0

    call system_clock(count_rate=cr)
    rate = real(cr)

   allocate(radint(natoms,nspecies,nnmax,lmax+1,nmax),&
     &     efact(natoms,nspecies,nnmax),length(natoms,nspecies,nnmax),nneigh(natoms,nspecies))

   nneigh(:,:) = 0
   length(:,:,:) = 0.d0
   efact(:,:,:) = 0.d0
   omega(:,:,:,:,:) = 0.d0
   harmonic(:,:,:,:,:) = 0.d0
   radint(:,:,:,:,:) = 0.d0
   orthoradint(:,:,:,:,:) = 0.d0

   call system_clock(ts)

   if (.not.periodic) then

    iat = 1
    ncentype = count(all_centres)
    ! Loop over species to centre on
    do icentype=1,nelements
     if (all_centres(icentype)) then
      ! Loop over centres of that species
      do icen=1,nneighmax(icentype)
       cen = all_indices(icentype,icen)
       ! Loop over all the species to use as neighbours
       k = 0
       do ispe=1,nelements
        if (all_species(ispe)) then
         k = k + 1
         ! Loop over neighbours of that species
         n = 1
         do ineigh=1,nneighmax(ispe)
          neigh = all_indices(ispe,ineigh)
          ! Compute distance vector
          rx = xyz(neigh,1) - xyz(cen,1)
          ry = xyz(neigh,2) - xyz(cen,2)
          rz = xyz(neigh,3) - xyz(cen,3)
          r2 = rx*rx + ry*ry + rz*rz
          ! Within cutoff?
          if (r2 .le. rcut2) then
           ! Central atom?
           if (neigh.eq.cen) then
            length(iat,k,n) = 0.d0
            efact(iat,k,n) = 1.d0
            harmonic(iat,k,0+1,0+1,n) = spherical_harmonic(0,0,0.d0,0.d0)
            nneigh(iat,k) = nneigh(iat,k) + 1
            n = n + 1
           else
            rdist = dsqrt(r2)
            length(iat,k,n) = rdist
            cth = rz/rdist
            ph = datan2(ry,rx)
            efact(iat,k,n) = dexp(-alpha*r2) * radial_scaling(rdist,radial_c,radial_r0,radial_m)
            do lval=0,lmax
             do im=0,2*lval
              mval = im-lval
              harmonic(iat,k,lval+1,im+1,n) = dconjg(spherical_harmonic(lval,mval,cth,ph))
             enddo
            enddo
            nneigh(iat,k) = nneigh(iat,k) + 1
            n = n + 1
           endif
          endif
         enddo
        endif
       enddo
       iat = iat + 1
      enddo
     endif
    enddo

   else

    ! Invert unit cell
    invcell(:,:) = cell(:,:)
    lwork = 1000
    call DGETRF(3,3,invcell,3,ipiv,info)
    call DGETRI(3,invcell,3,ipiv,work,lwork,info)

    iat = 1
    ncentype = count(all_centres)
    ! Loop over species to centre on
    do icentype=1,nelements
     if (all_centres(icentype)) then
      ! Loop over centres of that species
      do icen=1,nneighmax(icentype)
       cen = all_indices(icentype,icen)
       ! Loop over all the species to use as neighbours
       k = 0
       do ispe=1,nelements
        if (all_species(ispe)) then
         k = k + 1
         ! Loop over neighbours of that species
         n = 1
         do ineigh=1,nneighmax(ispe)
          neigh = all_indices(ispe,ineigh)
          ! Compute distance vector
          rv(1) = xyz(neigh,1) - xyz(cen,1)
          rv(2) = xyz(neigh,2) - xyz(cen,2)
          rv(3) = xyz(neigh,3) - xyz(cen,3)
          ! Apply periodic boundary conditions
          sv = matmul(invcell,rv)
          sv(1) = sv(1) - nint(sv(1))
          sv(2) = sv(2) - nint(sv(2))
          sv(3) = sv(3) - nint(sv(3))
          rcv = matmul(cell,sv)
          ! Loop over periodic images
          do ia=-ncell,ncell
           do ib=-ncell,ncell
            do ic=-ncell,ncell
             rx = rcv(1) + cell(1,1)*ia + cell(1,2)*ib + cell(1,3)*ic
             ry = rcv(2) + cell(2,1)*ia + cell(2,2)*ib + cell(2,3)*ic
             rz = rcv(3) + cell(3,1)*ia + cell(3,2)*ib + cell(3,3)*ic
             r2 = rx*rx + ry*ry + rz*rz
             ! Within cutoff?
             if (r2 .le. rcut2) then
              ! Central atom?
              if (neigh.eq.cen) then
               length(iat,k,n) = 0.d0
               efact(iat,k,n) = 1.d0
               harmonic(iat,k,0+1,0+1,n) = spherical_harmonic(0,0,0.d0,0.d0)
               nneigh(iat,k) = nneigh(iat,k) + 1
               n = n + 1
              else
               rdist = dsqrt(r2)
               length(iat,k,n) = rdist
               cth = rz/rdist
               ph = datan2(ry,rx)
               efact(iat,k,n) = dexp(-alpha*r2) * radial_scaling(rdist,radial_c,radial_r0,radial_m)
               do lval=0,lmax
                do im=0,2*lval
                 mval = im-lval
                 harmonic(iat,k,lval+1,im+1,n) = dconjg(spherical_harmonic(lval,mval,cth,ph))
                enddo
               enddo
               nneigh(iat,k) = nneigh(iat,k) + 1
               n = n + 1
              endif
             endif
            enddo
           enddo
          enddo
         enddo
        endif
       enddo
       iat = iat + 1
      enddo
     endif
    enddo

   endif

   call system_clock(tf)
!   write(*,'(A,F6.3,A)') ' got harmonic in ',(tf-ts)/rate,' s'
   call system_clock(ts)

   ! Get radial integral
   !$OMP PARALLEL DO SHARED(radint,efact,length,sg2) PRIVATE(n1,n2,normfact,sigmafact,i,k,nn,l,f1)
   do n1=1,nmax
    n2 = n1-1
    normfact = dsqrt(2.d0 / (gamma(1.5d0 + n2)*sigma(n1)**(3.d0 + 2.d0*n2)))
    sigmafact = (sg2**2 + sg2*sigma(n1)**2)/sigma(n1)**2
    do l=0,lmax
     f1 = 2.d0**(-0.5d0*(1.d0 + l-n2)) * (1.d0/sg2 + 1.d0/sigma(n1)**2)**(-0.5d0*(3.d0 + l+n2)) * &
     &   (gamma(0.5d0*(3.d0+l+n2))/gamma(1.5d0+l))  
     do nn=1,nnmax
      do k=1,nspecies
       do i=1,natoms
        radint(i,k,nn,l+1,n1) = efact(i,k,nn) * f1 * (length(i,k,nn)/sg2)**l * &
     &     hg(0.5d0*(3.d0+l+n2),1.5d0 + l,0.5d0 * length(i,k,nn)**2/sigmafact) 
        enddo
       enddo
     enddo
    enddo
    radint(:,:,:,:,n1) = radint(:,:,:,:,n1) * normfact
   enddo
   !$OMP END PARALLEL DO

   call system_clock(tf)
!   write(*,'(A,F6.3,A)') ' got radint in ',(tf-ts)/rate,' s'
   call system_clock(ts)

   ! Get orthoradint
   !$OMP PARALLEL DO SHARED(orthoradint,radint,orthomatrix) PRIVATE(iat,k,neigh)
   do l=1,lmax+1
    do k=1,nspecies
     do iat=1,natoms
      do neigh=1,nneigh(iat,k)
        orthoradint(iat,k,l,:,neigh) = matmul(orthomatrix,radint(iat,k,neigh,l,:))
       enddo
      enddo
    enddo
   enddo
   !$OMP END PARALLEL DO

   call system_clock(tf)
!   write(*,'(A,F6.3,A)') ' got orthoradint in ',(tf-ts)/rate,' s'
   call system_clock(ts)

   ! Get omega
   omega(:,:,:,:,:) = 0.d0
   !$OMP PARALLEL DO SHARED(omega,orthoradint,harmonic) PRIVATE(l,im,n1,k,iat,i)
   do im=1,2*lmax+1
     do n1=1,nmax
      do l=1,lmax+1
       do k=1,nspecies
        do iat=1,natoms
         omega(iat,k,n1,l,im) = dot_product(orthoradint(iat,k,l,n1,:),harmonic(iat,k,l,im,:))
        enddo
       enddo
      enddo
     enddo
   enddo
   !$OMP END PARALLEL DO

   call system_clock(tf)
!   write(*,'(A,F6.3,A)') ' got omega in ',(tf-ts)/rate,' s'

   deallocate(radint,efact,length,nneigh)

 end subroutine

!***************************************************************************************************

 real*8 function radial_scaling(r,c,r0,m)
  implicit none

   real*8 r,c,r0,m

   if (m.lt.1.d-10) then
    radial_scaling = 1.d0
    return
   endif
   if (c.eq.0.d0) then
    radial_scaling = 1.d0 / (r/r0)**m
   else
    radial_scaling = c / (c + (r/r0)**m)
   endif

 end function

!***************************************************************************************************

 function complex_to_real_matrix(degen)
  implicit none

   integer degen,lval,j
   complex*16 complex_to_real_matrix(degen,degen),st

    complex_to_real_matrix(:,:) = 0.d0
    lval = (degen-1)/2
    st = (-1.d0)**(lval+1)
    do j=0,lval-1
     complex_to_real_matrix(j+1,j+1) = (0.d0,1.d0)
     complex_to_real_matrix(j+1,degen-j) = st*(0.d0,1.d0)
     complex_to_real_matrix(degen-j,j+1) = 1.d0
     complex_to_real_matrix(degen-j,degen-j) = st*(-1.d0)
     st = st * (-1.d0)
    enddo
    complex_to_real_matrix(lval+1,lval+1) = dsqrt(2.d0)
    complex_to_real_matrix(:,:) = complex_to_real_matrix(:,:) / dsqrt(2.d0)
    complex_to_real_matrix(:,:) = dconjg(complex_to_real_matrix(:,:))

 end function

!***************************************************************************************************

! complex*16 function spherical_harmonic(l,m,costheta,phi)
!  implicit none
!
!   integer l,m,mm
!   real*8 costheta,phi
!   complex*16 rawfactor
!
!   mm = abs(m)
!   rawfactor = ((1.d0,0.d0)*cos(m*phi) + (0.d0,1.d0)*sin(m*phi)) * plgndr(l,abs(m),costheta)
!   spherical_harmonic = rawfactor * dsqrt( ((2*l + 1) / (4.d0*dacos(-1.d0)) * fact(l-mm)/fact(l+mm)))
!   if (m.lt.0) spherical_harmonic = spherical_harmonic * (-1.d0)**m
!
! end function

!***************************************************************************************************

! real*8 function plgndr(l,m,x)
!  implicit none
!
!   ! Subroutine from Numerical Recipes in Fortran
!
!   integer l,m
!   real*8 x
!   integer i,ll
!   real*8 fact,pll,pmm,pmmp1,somx2
!   if (m.lt.0.or.m.gt.l.or.abs(x).gt.1) stop 'ERROR: bad arguments in plgndr!'
!   pmm = 1.d0
!   if (m.gt.0) then
!    somx2 = sqrt((1.d0-x)*(1.d0+x))
!    fact = 1.d0
!    do i=1,m
!     pmm = -pmm*fact*somx2
!     fact = fact + 2.d0
!    enddo
!   endif
!   if (l.eq.m) then
!    plgndr=pmm
!   else
!    pmmp1 = x*(2*m+1)*pmm
!    if (l.eq.m+1) then
!     plgndr = pmmp1
!    else
!     do ll=m+2,l 
!      pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
!      pmm = pmmp1
!      pmmp1 = pll
!     enddo
!     plgndr = pll
!    endif
!   endif

!   return
! end function

!***************************************************************************************************

! real*8 function fact(n)
!  implicit none
!
!   integer n,i
!
!   if (n.lt.0) stop 'ERROR: positive number required for factorial!'
!   fact = 1.d0
!   do i=2,n
!    fact = fact * i
!   enddo
!
! end function

!***************************************************************************************************

 real*8 function maxdet(cell,nframes)
  implicit none

   integer nframes,i
   real*8 cell(nframes,3,3),cl(3,3),dt

   maxdet = 0.d0
   do i=1,nframes
    cl = cell(i,:,:)
    dt = cl(1,1)*cl(2,2)*cl(3,3) - cl(1,1)*cl(2,3)*cl(3,2) - cl(2,1)*cl(1,2)*cl(3,3) + &
     &     cl(1,2)*cl(2,3)*cl(3,1) + cl(1,3)*cl(2,1)*cl(3,2) - cl(1,3)*cl(2,2)*cl(3,1)
    maxdet = max(maxdet,dt)
   enddo

 end function

!***************************************************************************************************

end module
