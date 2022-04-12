module lode

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
   integer ineigh,neigh,lval,im,mval,n1,n2,l,i,k,nn,ncell,ia,ib,ic
   real*8 rcut2,rx,ry,rz,r2,rdist,cth,ph,normfact,sigmafact,rv(3),sv(3),rcv(3)
   complex*16 omega(natoms,nspecies,nmax,lmax+1,2*lmax+1)
   real*8 sg,rcut,xyz(natmax,3)
   real*8, allocatable :: radint(:,:,:,:,:),efact(:,:,:),length(:,:,:)
   real*8 sigma(nmax),orthomatrix(nmax,nmax),alpha,sg2,radial_c,radial_r0,radial_m,invcell(3,3)
   integer, allocatable :: nneigh(:,:)
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
!
!    # Define atomic grid for potential 
!    #radsize = 50 # number of radial points
!    gauss_points,gauss_weights = gausslegendre.gauss_legendre.gaulegf(x1=0.0,x2=rcut*2.0,n=radsize)
!   lebedev_grid = ld(lebsize)
!    spherical_grid = np.zeros((lebsize*radsize,3),float)
!    integration_weights = np.zeros((lebsize*radsize),float)
!    igrid = 0
!    for ir in range(radsize):
!        r = gauss_points[ir]
!        for ileb in range(lebsize):
!            spherical_grid[igrid] = r * lebedev_grid[ileb][:3]
!            integration_weights[igrid] = 4.0*np.pi * r**2 * lebedev_grid[ileb][3] * gauss_weights[ir] 
!            igrid += 1
!
!    # GET POLAR COORDINATES OVER ATOMIC GRID
!    lr = np.sqrt(np.sum(spherical_grid**2,axis=1))
!    lth = np.arccos(spherical_grid[:,2]/lr)
!    lph = np.arctan2(spherical_grid[:,1],spherical_grid[:,0])
!
!    harmonics = np.zeros(((lmax+1)**2,lebsize*radsize),complex) 
!    lm = 0
!    for l in range(lmax+1):
!        for im in range(2*l+1):
!            harmonics[lm,:] = np.conj(sc.sph_harm(im-l,l,lph[:],lth[:])) 
!            lm += 1
!
!    radial = radial_1D_mesh(sigma,nmax,lr,lebsize*radsize)
!    orthoradial = np.dot(orthomatrix,radial)
!
!    # compute near-field potential and project it onto the atomic basis
!    omega_near = nearfield.nearfield(nat,nspecies,nmax,lmax,lebsize*radsize,nneigh_near,alpha,coordx_near,spherical_grid,orthoradial,harmonics,integration_weights) 
!    omega_near = np.transpose(omega_near,(4,3,2,1,0))
!
!
!
!    return omega_near
!
!LIB = ctypes.CDLL(os.path.dirname(__file__) + "/lebedev.so")
!LDNS = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434,
!        590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890,
!        4334, 4802, 5294, 5810]
!
!def ld(n):
!    '''Returns (x, y, z) coordinates along with weight w. Choose among [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]'''
!
!    if n not in LDNS:
!        raise ValueError("n = {} not supported".format(n))
!    xyzw = np.zeros((4, n), dtype=np.float64)
!    c_double_p = ctypes.POINTER(ctypes.c_double)
!    n_out = ctypes.c_int(0)
!    getattr(LIB, "ld{:04}_".format(n))(
!        xyzw[0].ctypes.data_as(c_double_p),
!        xyzw[1].ctypes.data_as(c_double_p),
!        xyzw[2].ctypes.data_as(c_double_p),
!        xyzw[3].ctypes.data_as(c_double_p),
!        ctypes.byref(n_out))
!    assert n == n_out.value
!    return xyzw.T
!

 end subroutine

!***************************************************************************************************
 subroutine radial_1D_mesh(radial,sigma,nmax,rvec,rsize)
  implicit none

  ! Evaluate equispaced and normalied radial GTOs over a 1D mesh
  integer nmax,rsize,n,r
  integer rvec(rsize)
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
