module sagpr

 integer, parameter :: nelements = 118
 character(len=2), parameter :: atomic_names(nelements) = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si', &
     &     'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
     &     'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba', &
     &     'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
     &     'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
     &     'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)
 complex*16, allocatable :: PS(:,:,:,:)
 integer, allocatable :: components(:,:)
 real*8 start,finish

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
   k = 0
   do i=1,nmol
    do mu=1,degen
     k = k + 1
     do_prediction(i,mu) = prediction_lm(k)
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

 subroutine do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres, &
     &     all_species,ncut,sparsification,rs,periodic,mult_by_A)
  use SHTOOLS, only: wigner3j
  implicit none

  integer nframes,natmax,lm,nmax,lmax,ncut,degen,featsize,nnmax,nsmax,ispe,i,j,k,nspecies,n1,n2,l1,l2,l,mu,nu
  integer llmax,ps_shape(4),m,n,nn,jmin,jmax,im,ia,ib,mm
  real*8 xyz(nframes,natmax,3),rs(3),rcut,sg,cell(nframes,3,3)
  complex*16 sparsification(2,ncut,ncut)
  real*8 sigma(nmax),overlap(nmax,nmax),eigenval(nmax),diagsqrt(nmax,nmax),orthomatrix(nmax,nmax),inner
  character(len=4) atname(nframes,natmax)
  integer natoms(nframes),ipiv(nmax)
  logical periodic,all_centres(:),all_species(:),mult_by_A
  integer, allocatable :: all_indices(:,:,:),nneighmax(:,:),ncen(:),lvalues(:,:)
  integer, parameter :: lwmax = 10000
  integer info,lwork,work(lwmax)
  complex*16, allocatable :: omega(:,:,:,:,:),harmonic(:,:,:,:,:),omegatrue(:,:,:,:,:),omegaconj(:,:,:,:,:),ps_row(:,:,:)
  complex*16, allocatable :: CC(:,:),inner_mu(:,:)
  real*8, allocatable :: orthoradint(:,:,:,:,:),w3j(:,:,:,:),tmp_3j(:)
  integer, allocatable :: index_list(:)

  ! Get maximum number of neighbours
  if (.not.periodic) then
   nnmax = natmax
  else
   nnmax = natmax * 10
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
  if (.not.allocated(PS)) then
   allocate(PS(nframes,natmax,degen,featsize))
  else
   ! We have already allocated the power spectrum; let's see if we need to reallocate it
   ps_shape = shape(PS)
   if (.not.(ps_shape(1).ge.nframes .and. ps_shape(2).ge.natmax .and. ps_shape(3).ge.degen .and. ps_shape(4).ge.featsize)) then
    deallocate(PS)
    allocate(PS(nframes,natmax,degen,featsize))
   endif
  endif
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

  ! Do the power spectrum computation
  do i=1,nframes

   ! Get omega, harmonic and radint matrices
   allocate(omega(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),harmonic(natoms(i),nspecies,lmax+1,2*lmax+1,nnmax),&
     &     orthoradint(natoms(i),nspecies,lmax+1,nmax,nnmax))
   call initsoap(omega,harmonic,orthoradint,natoms(i),nspecies,nmax,lmax,nnmax,periodic,all_indices(i,:,:), &
     &     nneighmax(i,:),natmax,nsmax,cell(i,:,:),rs,sg,all_centres,all_species,rcut,xyz(i,:,:),sigma,orthomatrix)

   ! Compute power spectrum
   if (lm.eq.0) then
    ! Scalar
    allocate(omegatrue(natoms(i),nspecies,nmax,lmax+1,2*lmax+1),omegaconj(natoms(i),nspecies,nmax,lmax+1,2*lmax+1))
    do l=0,lmax
     omegatrue(:,:,:,l+1,:) = omega(:,:,:,l+1,:) / dsqrt(dsqrt(2.d0*l+1.d0))
    enddo
    omegaconj(:,:,:,:,:) = dconjg(omegatrue)
    if (ncut.gt.0) then
     do j=1,ncut
      do k=1,natoms(i)
       PS(i,k,1,j) = dot_product(omegatrue(k,components(j,1),components(j,3),components(j,5),:), &
     &     omegatrue(k,components(j,2),components(j,4),components(j,5),:))
      enddo
     enddo
    else
     stop 'ERROR: no sparsification information given; this is not recommended!'
    endif
    deallocate(omegatrue,omegaconj)

   else
     ! Spherical

     if (ncut.gt.0) then
      !$OMP PARALLEL DO SHARED(components,PS,omega,orthoradint,harmonic,w3j) PRIVATE(j,ia,ib,nn,mm,l1,l2,k,l,im,n)
      do j=1,ncut
       ia = components(j,1)
       ib = components(j,2)
       nn = components(j,3)
       mm = components(j,4)
       l1 = components(j,5)
       l2 = components(j,6)
       do k=0,2*lm
        do l=1,natoms(i)
         PS(i,l,k+1,j) = 0.d0
         do im=0,2*lmax
          if (abs(im-l1-k+lm).le.l2) then
           do n=1,nnmax
            PS(i,l,k+1,j) = PS(i,l,k+1,j) + &
     &          (omega(l,ia,nn,l1+1,im+1)*orthoradint(l,ib,l2+1,mm,n)* &
     &          conjg(harmonic(l,ib,l2+1,l2+im-l1-k+lm+1,n))*w3j(k+1,l1+1,l2+1,im+1))
           enddo
          endif
         enddo
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO
     else
      stop 'ERROR: no sparsification information given; this is not recommended!'
     endif
   endif

   ! Deallocate
   deallocate(omega,harmonic,orthoradint)

  enddo

  ! Multiply by A matrix
  if (mult_by_A) then
   do i=1,nframes
    do j=1,natmax
     do k=1,degen
      PS(i,j,k,:) = matmul(PS(i,j,k,:),sparsification(2,:,:))
     enddo
    enddo
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
  if (allocated(components)) deallocate(components)

 end subroutine

!***************************************************************************************************

 subroutine initsoap(omega,harmonic,orthoradint,natoms,nspecies,nmax,lmax,nnmax,periodic,all_indices,nneighmax, &
     &     natmax,nsmax,cell,rs,sg,all_centres,all_species,rcut,xyz,sigma,orthomatrix)
  implicit none

   integer natoms,nspecies,nmax,lmax,nnmax,natmax,nsmax,iat,ncentype,icentype,icen,cen,n,ispe
   integer ineigh,neigh,lval,im,mval,n1,n2,l,i,k,nn,ncell,ia,ib,ic
   real*8 rcut2,rx,ry,rz,r2,rdist,cth,ph,normfact,sigmafact,rv(3),sv(3),rcv(3)
   complex*16 omega(natoms,nspecies,nmax,lmax+1,2*lmax+1),harmonic(natoms,nspecies,lmax+1,2*lmax+1,nnmax)
   real*8 radint(natoms,nspecies,nnmax,lmax+1,nmax),orthoradint(natoms,nspecies,lmax+1,nmax,nnmax)
   real*8 efact(natoms,nspecies,nnmax),length(natoms,nspecies,nnmax),cell(3,3),rs(3),sg,rcut,xyz(natmax,3)
   real*8 sigma(nmax),orthomatrix(nmax,nmax),alpha,sg2,radial_c,radial_r0,radial_m,invcell(3,3)
   integer nneigh(natoms,nspecies),all_indices(nsmax,natmax),nneighmax(nsmax),ipiv(3),info,lwork,work(1000)
   logical periodic,all_centres(nelements),all_species(nelements)

   sg2 = sg*sg
   alpha = 1.d0 / (2.d0 * sg2)
   radial_c = rs(1)
   radial_r0 = rs(2)
   radial_m = rs(3)
   rcut2 = rcut*rcut
   ncell = 2

   nneigh(:,:) = 0
   length(:,:,:) = 0.d0
   efact(:,:,:) = 0.d0
   omega(:,:,:,:,:) = 0.d0
   harmonic(:,:,:,:,:) = 0.d0
   radint(:,:,:,:,:) = 0.d0
   orthoradint(:,:,:,:,:) = 0.d0

   call cpu_time(start)

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

   ! Get radial integral
   !$OMP PARALLEL DO SHARED(radint,efact,length,sg2) PRIVATE(n1,n2,normfact,sigmafact,i,k,nn,l)
   do n1=1,nmax
    n2 = n1-1
    normfact = dsqrt(2.d0 / (gamma(1.5d0 + n2)*sigma(n1)**(3.d0 + 2.d0*n2)))
    sigmafact = (sg2**2 + sg2*sigma(n1)**2)/sigma(n1)**2
    do l=0,lmax
     do nn=1,nnmax
      do k=1,nspecies
       do i=1,natoms
        radint(i,k,nn,l+1,n1) = efact(i,k,nn) * 2.d0**(-0.5d0*(1.d0 + l-n2)) * &
     &     (1.d0/sg2 + 1.d0/sigma(n1)**2)**(-0.5d0*(3.d0 + l+n2)) * &
     &     (gamma(0.5d0*(3.d0+l+n2))/gamma(1.5d0+l)) * (length(i,k,nn)/sg2)**l * &
     &     hg(0.5d0*(3.d0+l+n2),1.5d0 + l,0.5d0 * length(i,k,nn)**2/sigmafact) 
        enddo
       enddo
     enddo
    enddo
    radint(:,:,:,:,n1) = radint(:,:,:,:,n1) * normfact
   enddo
   !$OMP END PARALLEL DO

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

   ! Get omega
   omega(:,:,:,:,:) = 0.d0
   !$OMP PARALLEL DO SHARED(omega,orthoradint,harmonic) PRIVATE(l,im,n1,k,iat,i)
   do im=1,2*lmax+1
    do i=1,nnmax
     do n1=1,nmax
      do l=1,lmax+1
       do k=1,nspecies
        do iat=1,natoms
         omega(iat,k,n1,l,im) = omega(iat,k,n1,l,im) + (orthoradint(iat,k,l,n1,i) * harmonic(iat,k,l,im,i))
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
   !$OMP END PARALLEL DO

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

 complex*16 function spherical_harmonic(l,m,costheta,phi)
  implicit none

   integer l,m,mm
   real*8 costheta,phi
   complex*16 rawfactor

   mm = abs(m)
   rawfactor = ((1.d0,0.d0)*cos(m*phi) + (0.d0,1.d0)*sin(m*phi)) * plgndr(l,abs(m),costheta)
   spherical_harmonic = rawfactor * dsqrt( ((2*l + 1) / (4.d0*dacos(-1.d0)) * fact(l-mm)/fact(l+mm)))
   if (m.lt.0) spherical_harmonic = spherical_harmonic * (-1.d0)**m

 end function

!***************************************************************************************************

 real*8 function plgndr(l,m,x)
  implicit none

   ! Subroutine from Numerical Recipes in Fortran

   integer l,m
   real*8 x
   integer i,ll
   real*8 fact,pll,pmm,pmmp1,somx2
   if (m.lt.0.or.m.gt.l.or.abs(x).gt.1) stop 'ERROR: bad arguments in plgndr!'
   pmm = 1.d0
   if (m.gt.0) then
    somx2 = sqrt((1.d0-x)*(1.d0+x))
    fact = 1.d0
    do i=1,m
     pmm = -pmm*fact*somx2
     fact = fact + 2.d0
    enddo
   endif
   if (l.eq.m) then
    plgndr=pmm
   else
    pmmp1 = x*(2*m+1)*pmm
    if (l.eq.m+1) then
     plgndr = pmmp1
    else
     do ll=m+2,l 
      pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
      pmm = pmmp1
      pmmp1 = pll
     enddo
     plgndr = pll
    endif
   endif

   return
 end function

!***************************************************************************************************

 real*8 function fact(n)
  implicit none

   integer n,i

   if (n.lt.0) stop 'ERROR: positive number required for factorial!'
   fact = 1.d0
   do i=2,n
    fact = fact * i
   enddo

 end function

!***************************************************************************************************

 real*8 function hg ( a, b, x)

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) HG, the value of M(a,b,x).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) b
  real ( kind = 8 ) hg1
  real ( kind = 8 ) hg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rg
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) tba
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xg
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1

  pi = dacos(-1.d0)
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0D+00

  if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0D+300
  else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
    hg = 1.0D+00
  else if ( a == -1.0D+00 ) then
    hg = 1.0D+00 - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0D+00 ) then
    hg = ( 1.0D+00 + x / b ) * exp ( x )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
    hg = ( exp ( x ) - 1.0D+00 ) / x
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
    m = int ( - a )
    r = 1.0D+00
    hg = 1.0D+00
    do k = 1, m
      r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0D+00 ) then
    return
  end if

  if ( x < 0.0D+00 ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0D+00 ) then
    nl = 0
  end if

  if ( 2.0D+00 <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0D+00
  end if

  do n = 0, nl

    if ( 2.0D+00 <= a0 ) then
      a = a + 1.0D+00
    end if

    if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

      hg = 1.0D+00
      rg = 1.0D+00
      do j = 1, 500
        rg = rg * ( a + j - 1.0D+00 ) &
          / ( j * ( b + j - 1.0D+00 ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0D-15 ) then
          exit
        end if
      end do

    else

      ta = gamma(a)
      tb = gamma(b)
      xg = b - a
      tba = gamma(xg)
      sum1 = 1.0D+00
      sum2 = 1.0D+00
      r1 = 1.0D+00
      r2 = 1.0D+00
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0D+00 <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0D+00
    end do
  end if

  if ( x0 < 0.0D+00 ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
 end function


!***************************************************************************************************

end module
