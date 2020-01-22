module sagpr

 integer, parameter :: nelements = 86
 character(len=2), parameter :: atomic_names(nelements) = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si', &
     &     'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
     &     'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba', &
     &     'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
     &     'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'/)

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

 function do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres, &
     &     all_species,ncut,sparsification,rs,periodic)
  implicit none

  integer nframes,natmax,lm,nmax,lmax,ncut,degen,featsize,nnmax,nsmax,ispe,i,j,k,nspecies,n1,n2,l1,l2
  integer llmax
  real*8 xyz(nframes,natmax,3),rs(3),rcut,sg,cell(nframes,3,3)
  complex*16 sparsification(2,ncut,ncut)
  real*8 sigma(nmax),overlap(nmax,nmax),eigenval(nmax),diagsqrt(nmax,nmax),orthomatrix(nmax,nmax)
  character(len=4) atname(nframes,natmax)
  integer natoms(nframes),ipiv(nmax)
  logical periodic,all_centres(:),all_species(:)
  integer, allocatable :: all_indices(:,:,:),nneighmax(:,:),ncen(:),lvalues(:,:)
  integer, parameter :: lwmax = 10000
  integer info,lwork,work(lwmax)

  complex*16, allocatable :: PS(:,:,:,:)
  real*8, allocatable :: do_power_spectrum(:,:,:,:)

  ! Get maximum number of neighbours
  nnmax = natmax * 10

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
     do k=1,natmax
      !if (atn
     enddo
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
!  allocate(do_power_spectrum(nframes,natmax,degen,featsize))
  allocate(PS(nframes,natmax,degen,featsize))
  allocate(do_power_spectrum(nframes,natmax,degen,featsize))
  PS(:,:,:,:) = (0.d0,0.d0)
  do_power_spectrum(:,:,:,:) = 0.d0

  ! Do the power spectrum computation

  ! Normalize power spectrum

  deallocate(all_indices,nneighmax,ncen,PS)
  if (allocated(lvalues)) deallocate(lvalues)

 end function

!***************************************************************************************************


end module
