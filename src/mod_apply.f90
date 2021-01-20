module apply
 use sagpr

 ! Variables for getting coordinates
 real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
 character(len=4), allocatable :: atname(:,:)
 character(len=1000), allocatable :: comment(:)
 integer, allocatable :: natoms(:)
 integer nframes,nlines,natmax,nmol
 ! Variables for getting models
 integer nargs,ncut,nfeat,ncut0,nfeat0,reals,bytes
 real*8, allocatable :: PS_tr_lam(:,:,:,:),PS_tr_0(:,:,:,:),ker(:,:), ker_lm(:,:,:,:),natoms_tr(:)
 real*8, allocatable :: prediction_lm(:,:),wt(:)
 complex*16, allocatable :: sparsification(:,:,:),sparsification0(:,:,:)
 real*8 meanval
 logical do_scalar
 ! Parameters for PS and kernel building
 complex*16, allocatable :: PS(:,:,:,:),PS0(:,:,:,:)
 integer lm,nmax,lmax,zeta,degen
 real*8 rcut,sg,rs(3)
 logical periodic
 integer nmax0,lmax0
 real*8 rcut0,sg0,rs0(3)
 ! Defaults
 integer, parameter :: nmax_default = 8, lmax_default = 6
 real*8, parameter :: rcut_default = 4.d0,sg_default = 0.3d0,rs_default(3) = (/0.d0,0.d0,0.d0/)
 ! Committee models
 integer nw
 real*8, allocatable :: wt_c(:,:),meanval_c(:),prediction_lm_c(:,:,:)
 real*8 nu
 ! Atomic predictions
 logical atomic
 integer tot_natoms
 complex*16, allocatable :: PS_atomic(:,:,:,:),PS0_atomic(:,:,:,:)
 real*8, allocatable :: natoms_at(:)
 character(len=4), allocatable :: atname_at(:)
 ! Other
 logical verbose

 contains

!****************************************************************************************************************
subroutine set_defaults()
 implicit none

  lm = 0
  nmax = nmax_default
  lmax = lmax_default
  rcut = rcut_default
  sg = sg_default
  all_centres(:) = .false.
  all_species(:) = .false.
  rs = rs_default
  periodic = .false.
  zeta = 1

end subroutine

!****************************************************************************************************************
subroutine predict_frame(rate)
 implicit none

    integer ts,tf,i,k
    real*8 rate

    ! Get power spectrum
    if (.not.do_scalar) then
     call system_clock(ts)
     call do_power_spectrum(PS,xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,F6.3,A)') 'Got PS in',(tf-ts)/rate,' s'
    else
     call system_clock(ts)
     call do_power_spectrum(PS,xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',lm,' PS in',(tf-ts)/rate,' s'
     call system_clock(ts)
     call do_power_spectrum_scalar(PS0,xyz,atname,natoms,cell,nframes,natmax,0,nmax0,lmax0,rcut0,sg0, &
     &     ncut0,sparsification0,rs0,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',0,' PS in',(tf-ts)/rate,' s'
    endif

    ! If we have asked for atomic predictions, we need to convert to atomic
    ! power spectra
    if (atomic) then
     ! Here, rather than using the fact that we know what these quantities are,
     ! we will instead use the size function -- so as to make this agnostic of
     ! what we know (a good direction for future code)
     tot_natoms = 0
     do i=1,size(PS,1)
      tot_natoms = tot_natoms + natoms(i)
     enddo
     nframes = tot_natoms
     natmax = 1
     if (allocated(atname_at)) deallocate(atname_at)
     allocate(PS_atomic(tot_natoms,1,size(PS,3),size(PS,4)),natoms_at(tot_natoms),atname_at(tot_natoms))
     ! Populate atomic power spectrum array from original power spectrum
     k = 1
     do i=1,size(PS,1)
      PS_atomic(k:k+natoms(i),1,:,:) = PS(i,1:natoms(i),:,:)
      natoms_at(k:k+natoms(i)) = natoms(i)
      atname_at(k:k+natoms(i)) = atname(i,1:natoms(i))
      k = k + natoms(i)
     enddo
     ! Create new power spectrum array with the corrected shape
     deallocate(PS)
     allocate(PS(tot_natoms,1,size(PS_atomic,3),size(PS_atomic,4)))
     PS(:,:,:,:) = PS_atomic(:,:,:,:)
     deallocate(PS_atomic)
     if (do_scalar) then
      ! Repeat for scalar power spectrum if appropriate
      allocate(PS0_atomic(tot_natoms,1,size(PS0,3),size(PS0,4)))
      k = 1
      do i=1,size(PS0,1)
       PS0_atomic(k:k+natoms(i),1,:,:) = PS0(i,1:natoms(i),:,:)
       k = k + natoms(i)
      enddo
      deallocate(PS0)
      allocate(PS0(tot_natoms,1,size(PS0_atomic,3),size(PS0_atomic,4)))
      PS0(:,:,:,:) = PS0_atomic(:,:,:,:)
      deallocate(PS0_atomic)
     endif
     ! Create new number-of-atoms array with the corrected shape
     deallocate(natoms)
     allocate(natoms(tot_natoms))
     natoms(:) = 1
    endif

    ! Get kernel
    if (allocated(natoms_tr)) deallocate(natoms_tr)
    allocate(natoms_tr(nmol))
    natoms_tr(:) = 1.d0
    call system_clock(ts)
    if (.not.do_scalar) then
     if (lm.eq.0) then
      ker = do_scalar_kernel(real(PS),PS_tr_lam,nframes,nmol,nfeat,nfeat,natmax,1,dfloat(natoms),natoms_tr,zeta,.false.)
     else
      if (zeta.gt.1) stop 'ERROR: zeta>1 has been selected but no scalar power spectrum given!'
      ker_lm = do_linear_spherical_kernel(real(PS),PS_tr_lam,nframes,nmol,nfeat,nfeat, &
       &     natmax,1,dfloat(natoms),natoms_tr,zeta,.false.,degen)
     endif
    else
     if (lm.eq.0) stop 'ERROR: you have asked for a spherical power spectrum but provided lambda=0!'
     ker_lm = do_nonlinear_spherical_kernel(real(PS),PS_tr_lam,real(PS0),PS_tr_0,nframes,nmol,nfeat,nfeat, &
     &     nfeat0,nfeat0,natmax,1,dfloat(natoms),natoms_tr,zeta,.false.,degen)
    endif
    call system_clock(tf)
    if (verbose) write(*,'(A,F6.3,A)') 'Got kernel in ',(tf-ts)/rate,' s'

    ! Get predictions
    if (allocated(prediction_lm_c)) deallocate(prediction_lm_c)
    allocate(prediction_lm_c(nmol,degen,nw))
    if (lm.eq.0) then
     prediction_lm_c = do_prediction_c(ker,wt_c,meanval_c,degen,nframes,nmol,nw)
    else
     prediction_lm_c = do_prediction_c(ker_lm,wt_c,meanval_c,degen,nframes,nmol,nw)
    endif

    if (atomic) then
     ! Rescale predictions by the number of atoms, so that they are related to
     ! the frame predictions by a summation
     do i=1,size(prediction_lm_c,1)
      prediction_lm_c(i,:,:) = prediction_lm_c(i,:,:) / natoms_at(i)
     enddo
     deallocate(natoms_at)
    endif

end subroutine

!****************************************************************************************************************

end module
