module lode
 use gauss_legendre

 integer, parameter :: nelements = 118
 character(len=2), parameter :: atomic_names(nelements) = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si', &
     &     'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
     &     'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba', &
     &     'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
     &     'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
     &     'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

!***************************************************************************************************

 type LODE_Model
  ! Variables for running LODE calculations
  logical nonorm,fixed_cell
  real*8 sigewald
  integer radsize,lebsize
  ! G vectors for periodic LODE
  real*8, allocatable :: Gvec(:,:),Gval(:)
  real*8 Gcut,store_cell(3,3),icell(3,3),alphaewald
  integer nG
  real*8, allocatable :: orthoradint2(:,:,:)
  complex*16, allocatable :: harmonics2(:,:)
 end type LODE_Model

!***************************************************************************************************

 contains

!***************************************************************************************************

 subroutine direct_potential(omega,natoms,nspecies,nmax,lmax,nnmax,all_indices,nneighmax, &
     &     natmax,nsmax,sg,rcut,xyz,sigma,orthomatrix,all_species,all_centres,radsize,lebsize)
  implicit none

   integer natoms,nspecies,nmax,lmax,nnmax,natmax,nsmax,iat,ncentype,icentype,icen,cen,n,ispe
   integer ineigh,neigh,lval,im,mval,n1,n2,l,i,k,nn,ncell,ia,ib,ic,igrid,ir,ileb,lm,j,m,n_near
   real*8 rcut2,rx,ry,rz,r2,rdist,cth,ph,normfact,sigmafact,rv(3),sv(3),rcv(3)
   complex*16 omega(natoms,nspecies,nmax,lmax+1,2*lmax+1)
   complex*16 omega_near(2*lmax+1,lmax+1,nmax,nspecies,natoms)
   complex*16, allocatable :: harmonics(:,:)
   real*8 sg,rcut,xyz(natmax,3),r
   real*8, allocatable :: radint(:,:,:,:,:),efact(:,:,:),length(:,:,:),lebedev_grid(:,:),spherical_grid(:,:)
   real*8, allocatable :: integration_weights(:),gauss_points(:),gauss_weights(:),lr(:),lth(:),lph(:),radial(:,:)
   real*8, allocatable :: orthoradial(:,:),coordx_near(:,:,:,:)
   real*8 sigma(nmax),orthomatrix(nmax,nmax),alpha,sg2,radial_c,radial_r0,radial_m,invcell(3,3)
   integer, allocatable :: nneigh_near(:,:)
   integer all_indices(nsmax,natmax),nneighmax(nsmax),ipiv(3),info,lwork,work(1000)
   logical periodic
   integer ts,tf,cr,radsize,lebsize
   real*8 rate
   logical all_species(nelements),all_centres(nelements)

   alpha = 0.5d0 / (sg*sg)

   ! Process coordinates
   allocate(coordx_near(natoms,nspecies,natoms,3),nneigh_near(natoms,nspecies))
   coordx_near(:,:,:,:) = 0.d0
   nneigh_near(:,:) = 0.d0
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
!        n_near = 1
        do ineigh=1,nneighmax(ispe)
         neigh = all_indices(ispe,ineigh)
         coordx_near(iat,k,ineigh,:) = xyz(neigh,:) - xyz(cen,:)
!         n_near = n_near + 1
         nneigh_near(iat,k) = nneigh_near(iat,k) + 1
        enddo
       endif
      enddo
      iat = iat + 1
     enddo
    endif
   enddo

   ! Atomic grid for potential
   allocate(gauss_points(radsize),gauss_weights(radsize))
   call gaulegf(0.d0,2.d0*rcut,gauss_points,gauss_weights,radsize)

   ! Get Lebedev grid
   allocate(lebedev_grid(4,lebsize))
   lebedev_grid(:,:) = 0.d0
   call get_lebedev_grid(lebedev_grid,lebsize)

   ! Get spherical grid and integration weights
   allocate(spherical_grid(lebsize*radsize,3),integration_weights(lebsize*radsize))
   spherical_grid(:,:) = 0.d0
   integration_weights(:) = 0.d0
   igrid = 1
   do ir=1,radsize
    r = gauss_points(ir)
    do ileb=1,lebsize
     spherical_grid(igrid,:) = r * lebedev_grid(1:3,ileb)
     integration_weights(igrid) = 4.d0 * dacos(-1.d0) * r * r * lebedev_grid(4,ileb) * gauss_weights(ir)
     igrid = igrid + 1
    enddo
   enddo

   ! Get polar coordinates on atomic grid
   allocate(lr(lebsize*radsize),lth(lebsize*radsize),lph(lebsize*radsize))
   lr(:) = 0.d0
   do i=1,lebsize*radsize
    lr(i) = dsqrt(dot_product(spherical_grid(i,:),spherical_grid(i,:)))
    lth(i) = dacos(spherical_grid(i,3)/lr(i))
    lph(i) = datan2(spherical_grid(i,2),spherical_grid(i,1))
   enddo

   ! Get spherical harmonics array
   allocate(harmonics((lmax+1)*(lmax+1),lebsize*radsize))
   harmonics(:,:) = (0.d0,0.d0)
   lm = 1
   do l=0,lmax
    do im=-l,l
     do i=1,lebsize*radsize
      harmonics(lm,i) = dconjg(spherical_harmonic(l,im,dcos(lth(i)),lph(i)))
     enddo
     lm = lm + 1
    enddo
   enddo

   ! Get orthoradial array
   allocate(radial(nmax,lebsize*radsize))
   call radial_1D_mesh(radial,sigma,nmax,lr,lebsize*radsize)
   allocate(orthoradial(nmax,lebsize*radsize))
   orthoradial = matmul(orthomatrix,radial)

   ! Compute near-field potential and project it onto the atomic basis
   call nearfield(natoms,nspecies,nmax,lmax,lebsize*radsize,nneigh_near,alpha,coordx_near,spherical_grid,orthoradial, &
     &     harmonics,integration_weights,omega)

   deallocate(lebedev_grid,spherical_grid,gauss_points,gauss_weights,lr,lth,lph,harmonics,radial,orthoradial)
   deallocate(coordx_near,nneigh_near)

 end subroutine
!***************************************************************************************************
 subroutine get_Gvectors(this,cell,invcell,nmax,lmax,rc,sigma,orthomatrix)
  implicit none

   type(LODE_Model) :: this
   integer nmax,lmax,i,j,k
   real*8 cell(3,3),invcell(3,3),radint(lmax+1,nmax),prefacts(nmax,lmax+1),M(3,3),detM
   real*8 rc,sigma(nmax),orthomatrix(nmax,nmax),normcell(3,3),Gvec_new(3)
   real*8, allocatable :: G_intermediate(:,:)
   integer n1max,n2max,n3max,numtot,n1,n2,n3

   ! First, check whether we already have G-vectors allocated
   if (allocated(this%Gvec)) then
    ! If the G-vectors have already been allocated, check whether the cell is the same as the previous one
    write(*,*) 'G VECTORS ALREADY ALLOCATED'
    if (maxval(abs(cell-this%store_cell)).gt.1.d-8) then
     write(*,*) 'G VECTORS MUST BE RECALCULATED'
    else
     write(*,*) 'G VECTORS ARE OK'
     return
    endif
   endif

   write(*,*) 'HAVE TO ALLOCATE G VECTORS'
   ! Get G vectors
   this%Gcut = 2.d0 * dacos(-1.d0) / (2.d0 * this%sigewald)
   call old_Gvec_generator(invcell,this%Gcut,this%nG,this%Gvec,this%Gval)

   ! Compute Fourier integrals for the cell
   this%alphaewald = 0.5d0 / (this%sigewald*this%sigewald)
   allocate(this%orthoradint2(lmax+1,nmax,this%nG),this%harmonics2(this%nG,(lmax+1)*(lmax+1))) 
   call fourier_integrals(this%nG,nmax,lmax,this%alphaewald,rc,sigma,this%Gval,this%Gvec, &
     &     radint,orthomatrix,prefacts,this%orthoradint2,this%harmonics2)

 end subroutine
!***************************************************************************************************
 subroutine old_Gvec_generator(invcell,Gcut,nG,Gvec,Gval)
  implicit none

   ! Gvec generator following the version in TENSOAP
   real*8 invcell(3,3),normcell(3,3),M(3,3),detM,Gcut,Gvec_new(3),Gcut2
   integer n1max,n2max,n3max,numtot,nG,n1,n2,n3,i
   real*8, allocatable :: G_intermediate(:,:),Gvec(:,:),Gval(:)

   ! Outer product matrix of cell vectors
   normcell = invcell * 2.d0 * dacos(-1.d0)
   M = matmul(normcell,transpose(normcell))
   ! Get determinant of matrix
   detM = M(1,1)*M(2,2)*M(3,3) + 2.d0*M(1,2)*M(2,3)*M(3,1) &
     &     -M(1,1)*M(2,3)*M(2,3) - M(2,2)*M(1,3)*M(1,3) &
     &     -M(3,3)*M(1,2)*M(1,2)
   ! Maximum numbers in each direction
   n1max = floor(dsqrt((M(2,2)*M(3,3)-M(2,3)*M(2,3))/detM)*Gcut)
   n2max = floor(dsqrt((M(1,1)*M(3,3)-M(1,3)*M(1,3))/detM)*Gcut)
   n3max = floor(dsqrt((M(1,1)*M(2,2)-M(1,2)*M(1,2))/detM)*Gcut)
   ! Estimate for total number of G vectors
   numtot = 1 + n3max + n2max*(2*n3max+1) + n1max*(2*n2max+1)*(2*n3max+1)
   allocate(G_intermediate(numtot,3))
   G_intermediate(:,:) = 0.d0

   ! Fill intermediate set of G vectors
   nG = 1
   G_intermediate(nG,:) = (/0.d0,0.d0,0.d0/)
   Gcut2 = Gcut*Gcut
   ! Step 1: points of the form (0,0,n3>0)
   do n3=1,n3max
    Gvec_new(:) = n3 * normcell(3,:)
    if (Gvec_new(1)*Gvec_new(1) + Gvec_new(2)*Gvec_new(2) + Gvec_new(3)*Gvec_new(3) .lt. Gcut*Gcut) then
     nG = ng + 1
     G_intermediate(nG,:) = Gvec_new(:)
    endif
   enddo
   ! Step 2: points of the form (0,n2>0,n3)
   do n2=1,n2max
    do n3=-n3max,n3max
     Gvec_new(:) = n2 * normcell(2,:) + n3*normcell(3,:)
     if (Gvec_new(1)*Gvec_new(1) + Gvec_new(2)*Gvec_new(2) + Gvec_new(3)*Gvec_new(3) .lt. Gcut*Gcut) then
      nG = nG + 1
      G_intermediate(nG,:) = Gvec_new(:)
     endif
    enddo
   enddo
   ! Step 3: remaining points of the form (n1>0,n2,n3)
   do n1=1,n1max
    do n2=-n2max,n2max
     do n3=-n3max,n3max
      Gvec_new(:) = n1*normcell(1,:) + n2*normcell(2,:) + n3*normcell(3,:)
      if (Gvec_new(1)*Gvec_new(1) + Gvec_new(2)*Gvec_new(2) + Gvec_new(3)*Gvec_new(3) .lt. Gcut*Gcut) then
       nG = nG + 1
       G_intermediate(nG,:) = Gvec_new(:)
      endif
     enddo
    enddo
   enddo
   ! Take only the nonzero points from this array to give Gvec
   ! Gval contains norms
   allocate(Gvec(nG,3),Gval(nG))
   do i=1,nG
    Gvec(i,:) = G_intermediate(i,:)
    Gval(i) = norm2(Gvec(i,:))
   enddo

   deallocate(G_intermediate)

 end subroutine
!***************************************************************************************************
 subroutine new_Gvec_generator(invcell,Gcut,nG,Gvec,Gval)
  implicit none

   ! Gvec generator following the C code by Andrea Grisafi and Kevin Kazuki Huguenin-Dumittan
   real*8 invcell(3,3),normcell(3,3),M(3,3),detM,Gcut,Gvec_new(3)
   integer n1max,n2max,n3max,numtot,nG,n1,n2,n3,i
   real*8, allocatable :: G_intermediate(:,:),Gvec(:,:),Gval(:)

   ! Outer product matrix of cell vectors
   normcell = invcell * 2.d0 * dacos(-1.d0)
   M = matmul(normcell,transpose(normcell))
   ! Get determinant of matrix
   detM = M(1,1)*M(2,2)*M(3,3) + 2.d0*M(1,2)*M(2,3)*M(3,1) &
     &     -M(1,1)*M(2,3)*M(2,3) - M(2,2)*M(1,3)*M(1,3) &
     &     -M(3,3)*M(1,2)*M(1,2)
   ! Maximum numbers in each direction
   n1max = floor(dsqrt((M(2,2)*M(3,3)-M(2,3)*M(2,3))/detM)*Gcut)
   n2max = floor(dsqrt((M(1,1)*M(3,3)-M(1,3)*M(1,3))/detM)*Gcut)
   n3max = floor(dsqrt((M(1,1)*M(2,2)-M(1,2)*M(1,2))/detM)*Gcut)
   ! Estimate for total number of G vectors
   numtot = 1 + n3max + n2max*(2*n3max+1) + n1max*(2*n2max+1)*(2*n3max+1)
   allocate(G_intermediate(numtot,3))
   G_intermediate(:,:) = 0.d0

   ! Fill intermediate set of G vectors
   nG = 1
   G_intermediate(nG,:) = (/0.d0,0.d0,0.d0/)
   Gvec_new(:) = (/0.d0,0.d0,0.d0/)
   ! Step 1: points of the form (0,0,n3>0)
   do n3=1,n3max
    Gvec_new(:) = Gvec_new(:) + normcell(3,:)
    if (norm2(Gvec_new).le.Gcut) then
     nG = ng + 1
     G_intermediate(nG,:) = Gvec_new(:)
    endif
   enddo
   ! Step 2: points of the form (0,n2>0,n3)
   do n2=1,n2max
    ! Update current vector for new n2 value.
    ! We subtract (n3max+1)*normcell(3,:) so that we only have to add normcell(3,:)
    ! at each iteration to get the correct vector
    Gvec_new(:) = n2*normcell(2,:) - (n3max+1)*normcell(3,:)
    do n3=0,2*n3max+1
     Gvec_new(:) = Gvec_new(:) + normcell(3,:)
     if (norm2(Gvec_new).le.Gcut) then
      nG = nG + 1
      G_intermediate(nG,:) = Gvec_new(:)
     endif
    enddo
   enddo
   ! Step 3: remaining points of the form (n1>0,n2,n3)
   do n1=1,n1max
    do n2=0,2*n2max+1
     ! Update current vector for new n2 value
     ! We subtract (n3max+1)*normcell(3,:) so that we only have to add normcell(3,:)
     ! at each iteration to get the correct vector
     Gvec_new(:) = n1*normcell(1,:) + n2*normcell(2,:) - n2max*normcell(2,:) - (n3max+1)*normcell(3,:)
     do n3=0,2*n3max+1
      Gvec_new(:) = Gvec_new(:) + normcell(3,:)
      if (norm2(Gvec_new).le.Gcut) then
       nG = nG + 1
       G_intermediate(nG,:) = Gvec_new(:)
      endif
     enddo
    enddo
   enddo
   ! Take only the nonzero points from this array to give Gvec
   ! Gval contains norms
   allocate(Gvec(nG,3),Gval(nG))
   do i=1,nG
    Gvec(i,:) = G_intermediate(i,:)
    Gval(i) = norm2(Gvec(i,:))
   enddo

   deallocate(G_intermediate)

 end subroutine
!***************************************************************************************************
 subroutine fourier_integrals(nG,nmax,lmax,alpha,rc,sigma,Gval,Gvec,radint,orthomatrix,prefacts,orthoradint,harmonics)
  implicit none

   integer nG,nmax,lmax,n,l,iG,lm,n1,n2,im
   real*8 radint(lmax+1,nmax),prefacts(nmax,lmax+1),normfact,G2,fourierpot
   real*8 alpha,rc,sigma(nmax),orthomatrix(nmax,nmax),th,ph,arghyper
   real*8 orthoradint(lmax+1,nmax,nG),Gval(nG),Gvec(nG,3)
   complex*16 harmonics(nG,(lmax+1)*(lmax+1))

   orthoradint(:,:,:) = 0.d0
   harmonics(:,:) = (0.d0,0.d0)

   ! Normalization factor for primitive radial functions
   do n=0,nmax-1
    normfact = dsqrt(2.d0 / (gamma(1.5d0 + n)*sigma(n+1)**(3.d0 + 2.d0*n)))
    do l=0,lmax
     ! Precompute common prefactors for each (n,l) pair
     prefacts(n+1,l+1) = normfact * dsqrt(dacos(-1.d0)) * 2.d0**( (n-l-1)/2.d0) * sigma(n+1)**(3.d0 + l + n) &
     &      * gamma(0.5d0*(3.d0+l+n)) / gamma(1.5d0 + l)
    enddo
   enddo

   do iG=2,nG
    G2 = Gval(iG)*Gval(iG)
    fourierpot = dexp(-G2 / (4.d0*alpha)) / G2
    do n=0,nmax-1
     arghyper = -0.5d0 * G2 * (sigma(n+1)*sigma(n+1))
     do l=0,lmax
      ! Radial integral
      radint(l+1,n+1) = prefacts(n+1,l+1) * fourierpot * Gval(iG)**l * &
     &     hg(0.5d0*(3.d0+l+n),1.5d0+l,arghyper)
     enddo
    enddo
    ! Compute polar angles
    th = dacos(Gvec(iG,3)/Gval(iG))
    ph = datan2(Gvec(iG,2),Gvec(iG,1))
    lm = 1
    do l=0,lmax
     do n1=0,nmax-1
      do n2=0,nmax-1
       ! Orthogonalize radial integrals with Loewdin
       orthoradint(l+1,n1+1,iG) = orthoradint(l+1,n1+1,iG) + orthomatrix(n1+1,n2+1)*radint(l+1,n2+1)
      enddo
     enddo
     do im=-l,l
       harmonics(iG,lm) = dconjg(spherical_harmonic(l,im,dcos(th),ph)) * (0.d0,1.d0)**l
       lm = lm + 1
     enddo
    enddo
   enddo

 end subroutine
!***************************************************************************************************
 subroutine ewald_potential(omega2,natoms,nspecies,nmax,lmax,nnmax,all_indices, &
     &     nneighmax,all_species,all_centres,rs,sg,rcut,xyz,cell,orthomatrix, &
     &     radsize,lebsize,sigewald,Gvec,Gval,nG,orthoradint,harmonics,sigma,natmax,nsmax)
  implicit none

   integer natoms,nmax,lmax,radsize,lebsize,nG,nspecies,natmax,nnmax,nsmax
   integer all_indices(nsmax,natmax),nneighmax(nsmax)
   real*8 Gvec(nG,3),Gval(nG),cell(3,3),xyz(natmax,3),sigewald
   real*8 sigma(nmax),orthomatrix(nmax,nmax),rs(3)
   complex*16 omega2(natoms,nspecies,nmax,lmax+1,2*lmax+1)
   logical all_species(nelements),all_centres(nelements)
   real*8 orthoradint(lmax+1,nmax,nG),sg,rcut
   complex*16 harmonics(nG,(lmax+1)*(lmax+1))

 end subroutine
!***************************************************************************************************
 subroutine get_lebedev_grid(lebedev_grid,grid_size)
  implicit none

  ! Get Lebedev integration grid
  integer grid_size
  real*8 lebedev_grid(4,grid_size)

  select case (grid_size)
   case (6)
    call LD0006(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (14)
    call LD0014(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (26)
    call LD0026(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (38)
    call LD0038(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (50)
    call LD0050(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (74)
    call LD0074(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (86) 
    call LD0086(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size) 
   case (110)
    call LD0110(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (146)
    call LD0146(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (170)
    call LD0170(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (194)
    call LD0194(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (230)
    call LD0230(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (266)
    call LD0266(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (302)
    call LD0302(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (350)
    call LD0350(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (434)
    call LD0434(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (590)
    call LD0590(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (770)
    call LD0770(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (974)
    call LD0974(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (1202)
    call LD1202(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (1454)
    call LD1454(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (1730)
    call LD1730(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (2030)
    call LD2030(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (2354)
    call LD2354(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (2702)
    call LD2702(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (3074)
    call LD3074(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (3470)
    call LD3470(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (3890)
    call LD3890(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (4334)
    call LD4334(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (4802)
    call LD4802(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (5294)
    call LD5294(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case (5810)
    call LD5810(lebedev_grid(1,:),lebedev_grid(2,:),lebedev_grid(3,:),lebedev_grid(4,:),grid_size)
   case default
    write(*,*) 'ERROR: incorrect Lebedev grid size!',grid_size
    stop
  end select

 end subroutine

!***************************************************************************************************
 subroutine radial_1D_mesh(radial,sigma,nmax,rvec,rsize)
  implicit none

  ! Evaluate equispaced and normalied radial GTOs over a 1D mesh
  integer nmax,rsize,n,r
  real*8 rvec(rsize)
  real*8 sigma(nmax),radial(nmax,rsize),inner,arg,rval

  do n=0,nmax-1
   inner = 0.5d0 * gamma(n + 1.5d0) * (sigma(n+1)*sigma(n+1))**(n + 1.5d0)
   do r=1,rsize
    rval = rvec(r)
    arg = -0.5d0 * rval*rval / (sigma(n+1)*sigma(n+1))
    radial(n+1,r) = (rval**n) * dexp(arg) / dsqrt(inner)
   enddo
  enddo
  

 end subroutine
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
