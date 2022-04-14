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
  real*8, allocatable :: Gvec(:,:),Gval(:),G(:),Gx(:,:)
  real*8 Gcut,store_cell(3,3),icell(3,3),alphaewald
  integer nside(3),nG,irad(3)
  integer, allocatable :: iGvec(:,:),imGvec(:,:),iGx(:,:),iGmx(:,:)
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
   integer nmax,lmax,i,j
   real*8 cell(3,3),invcell(3,3),radint(lmax+1,nmax),prefacts(nmax,lmax+1)
   real*8 rc,sigma(nmax),orthomatrix(nmax,nmax)

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

!   ! Wave-vectors for reciprocal space calculation
!   this%Gcut = 2.d0 * dacos(-1.d0) / (2.d0 * this%sigewald)
!   ! Allocate default for visualization grid
!   this%nside = (/256,256,256/)
!   ! Get inverse cell
!   do i=1,3
!    do j=1,3
!     this%icell(i,j) = invcell(i,j) * 2.d0 * dacos(-1.d0)
!    enddo
!   enddo

   ! Get wave-vectors
!   allocate(this%G(10000000),this%Gx(3,10000000),this%iGx(3,10000000))
!   allocate(this%iGmx(3,10000000))
!   call ggen(this%icell(1,:),this%icell(2,:),this%icell(3,:),this%nside,this%Gcut, &
!     & this%G,this%Gx,this%iGx,this%iGmx,this%nG,this%irad)

!   if ((this%nside(1)/2 .le. this%irad(1)) .or. (this%nside(2)/2 .le. this%irad(2)) &
!     &     .or. (this%nside(3)/2 .le. this%irad(3))) then
!    stop 'ERROR: G-ellipsoid covers more points than half of the box side: decrease Gcut or &
!     &use more grid points for 3D visualization'
!   endif

!   ! G moduli of  the semi-sphere x>0
!   allocate(this%Gval(this%nG))
!   this%Gval(:) = this%G(1:this%nG)
!   ! G vectors of the semi-sphere x>0
!   allocate(this%Gvec(this%nG,3))
!   do i=1,this%nG
!    do j=1,3
!     this%Gvec(i,j) = this%Gx(j,i)
!    enddo
!   enddo
!   ! G vector indices of the semi-sphere x>0
!   allocate(this%iGvec(this%nG,3))
!   do i=1,this%nG
!    do j=1,3
!     this%iGvec(i,j) = this%iGx(j,i)
!    enddo
!   enddo
!   ! G vector indices of the semi-sphere x<0
!   allocate(this%imGvec(this%nG,3))
!   do i=1,this%nG
!    do j=1,3
!     this%imGvec(i,j) = this%iGmx(j,i)
!    enddo
!   enddo

   ! Compute Fourier integrals for the cell
!   this%alphaewald = 0.5d0 / (this%sigewald*this%sigewald)
!   allocate(this%orthoradint(lmax+1,nmax,this%nG),this%harmonics(this%nG,(lmax+1)*(lmax+1))) 
!   call fourier_integrals(this%nG,nmax,lmax,this%alphaewald,rc,sigma,this%Gval,this%Gvec, &
!     &     radint,orthomatrix,prefacts,this%orthoradint,this%harmonics)

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

   do iG=1,nG
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
