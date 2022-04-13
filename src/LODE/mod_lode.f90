module lode
 use gauss_legendre

 integer, parameter :: nelements = 118
 character(len=2), parameter :: atomic_names(nelements) = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si', &
     &     'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
     &     'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba', &
     &     'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
     &     'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
     &     'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

 type LODE_Model
  ! Variables for running LODE calculations
  logical nonorm,fixed_cell
  real*8 sigewald
  integer radsize,lebsize
 end type LODE_Model

 contains

!***************************************************************************************************

 subroutine direct_potential(omega,natoms,nspecies,nmax,lmax,nnmax,all_indices,nneighmax, &
     &     natmax,nsmax,sg,rcut,xyz,sigma,orthomatrix,all_species,all_centres,radsize,lebsize)
  implicit none

   integer natoms,nspecies,nmax,lmax,nnmax,natmax,nsmax,iat,ncentype,icentype,icen,cen,n,ispe
   integer ineigh,neigh,lval,im,mval,n1,n2,l,i,k,nn,ncell,ia,ib,ic,igrid,ir,ileb,lm,j,m
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
!
!    species = np.zeros(nat,float)
!
!    # process coordinates 
!    coordx_near = np.zeros((nat,nspecies,nat,3), dtype=float)
!    nneigh_near = np.zeros((nat,nspecies),int)
!    iat = 0
!    ncentype = len(centers)
!    # loop over species to center on
!    for icentype in range(ncentype):
!        centype = centers[icentype]
!        # loop over centers of that species
!        for icen in range(nneighmax[centype]):
!            cen = atom_indexes[centype,icen]
!            # loop over all the species to use as neighbours
!            for ispe in range(nspecies):
!                spe = all_species[ispe]
!                # loop over neighbours of that species
!                n_near = 0
!                for ineigh in range(nneighmax[spe]):
!                    neigh = atom_indexes[spe,ineigh]
!                    rx = coords[neigh,0] - coords[cen,0]
!                    ry = coords[neigh,1] - coords[cen,1]
!                    rz = coords[neigh,2] - coords[cen,2]
!                    coordx_near[iat,ispe,n_near,0] = rx
!                    coordx_near[iat,ispe,n_near,1] = ry
!                    coordx_near[iat,ispe,n_near,2] = rz
!                    n_near += 1
!                    nneigh_near[iat,ispe] += 1
!            species[iat] = centype
!            iat = iat + 1

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
     spherical_grid(igrid,:) = r * lebedev_grid(ileb,:)
     integration_weights(igrid) = 4.d0 * dacos(-1.d0) * r * r * lebedev_grid(ileb,3) * gauss_weights(ir)
     igrid = igrid + 1
    enddo
   enddo

   ! Get polar coordinates on atomic grid
   allocate(lr(lebsize*radsize),lth(lebsize*radsize),lph(lebsize*radsize))
   lr(:) = 0.d0
   do i=1,lebsize*radsize
    lr(i) = dot_product(spherical_grid(i,:),spherical_grid(i,:))
    lth(i) = dacos(spherical_grid(i,3)/lr(i))
    lph(i) = datan2(spherical_grid(i,2),spherical_grid(i,1))
   enddo

   ! Get spherical harmonics array
   allocate(harmonics((lmax+1)*(lmax+1),lebsize*radsize))
   harmonics(:,:) = (0.d0,0.d0)
   lm = 1
   do l=0,lmax
    do im=1,2*l+1
     do i=1,lebsize*radsize
      harmonics(lm,i) = dconjg(spherical_harmonic(l,im-l,dcos(lth(i)),lph(i)))
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
     &     harmonics,integration_weights,omega_near)
   do i=1,natoms
    do j=1,nspecies
     do k=1,nmax
      do l=1,lmax+1
       do m=1,2*lmax+1
        omega(i,j,k,l,m) = omega_near(m,l,k,j,i)
       enddo
      enddo
     enddo
    enddo
   enddo

   deallocate(lebedev_grid,spherical_grid,gauss_points,gauss_weights,lr,lth,lph,harmonics,radial,orthoradial)

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

end module
