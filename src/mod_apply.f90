module apply
 use sagpr

  ! Defaults
  integer, parameter :: nmax_default = 8, lmax_default = 6
  real*8, parameter :: rcut_default = 4.d0,sg_default = 0.3d0,rs_default(3) = (/0.d0,0.d0,0.d0/)

 type Frame_XYZ
  ! Variables for geting coordinates
  real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
  character(len=4), allocatable :: atname(:,:)
  character(len=1000), allocatable :: comment(:)
  integer, allocatable :: natoms(:)
  integer nframes,nlines,natmax
 end type Frame_XYZ

 type SAGPR_Model
  ! Variables for getting models
  integer ncut,nfeat,ncut0,nfeat0,nmol
  real*8, allocatable :: PS_tr_lam(:,:,:,:),PS_tr_0(:,:,:,:),ker(:,:), ker_lm(:,:,:,:),natoms_tr(:)
  real*8, allocatable :: prediction_lm(:,:),wt(:)
  complex*16, allocatable :: sparsification(:,:,:),sparsification0(:,:,:)
  real*8 meanval
  logical do_scalar
  ! Parameters for PS and kernel building
  complex*16, allocatable :: PS(:,:,:,:),PS0(:,:,:,:)
  integer, allocatable :: components(:,:),components0(:,:)
  real*8, allocatable :: w3j(:,:,:,:)
  integer lm,nmax,lmax,zeta,degen
  real*8 rcut,sg,rs(3)
  logical periodic
  integer nmax0,lmax0
  real*8 rcut0,sg0,rs0(3)
  logical all_species(nelements),all_centres(nelements)
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
 end type SAGPR_Model

 contains

!****************************************************************************************************************
subroutine set_defaults(this)
 implicit none

  type(SAGPR_Model), intent(inout) :: this

  this%lm = 0
  this%nmax = nmax_default
  this%lmax = lmax_default
  this%rcut = rcut_default
  this%sg = sg_default
  this%all_centres(:) = .false.
  this%all_species(:) = .false.
  this%rs = rs_default
  this%periodic = .false.
  this%zeta = 1

end subroutine
!****************************************************************************************************************
subroutine predict_frame(this,frames,rate)
 implicit none

    type(SAGPR_Model), intent(inout) :: this
    type(Frame_XYZ), intent(inout) :: frames

    integer ts,tf,i,k,backup_nframes,backup_natmax
    integer, allocatable :: natoms_backup(:)
    real*8 rate

    ! Get power spectrum
    if (.not.this%do_scalar) then
     call system_clock(ts)
     call do_power_spectrum(this%PS,this%components,frames%xyz,frames%atname, &
     &     frames%natoms,frames%cell,frames%nframes,frames%natmax,this%lm,this%nmax, &
     &     this%lmax,this%rcut,this%sg, &
     &     this%ncut,this%sparsification,this%rs,this%periodic,.true., &
     &     this%w3j,this%all_species,this%all_centres)
     call system_clock(tf)
     if (this%verbose) write(*,'(A,F6.3,A)') 'Got PS in',(tf-ts)/rate,' s'
    else
     call system_clock(ts)
     call do_power_spectrum(this%PS,this%components,frames%xyz,frames%atname, &
     &     frames%natoms,frames%cell,frames%nframes,frames%natmax,this%lm,this%nmax, &
     &     this%lmax,this%rcut,this%sg, &
     &     this%ncut,this%sparsification,this%rs,this%periodic,.true., &
     &     this%w3j,this%all_species,this%all_centres)
     call system_clock(tf)
     if (this%verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',this%lm,' PS in',(tf-ts)/rate,' s'
     call system_clock(ts)
     call do_power_spectrum_scalar(this%PS0,this%components0,frames%xyz,frames%atname, &
     &     frames%natoms,frames%cell,frames%nframes,frames%natmax,0,this%nmax0, &
     &     this%lmax0,this%rcut0,this%sg0, &
     &     this%ncut0,this%sparsification0,this%rs0,this%periodic,.true., &
     &     this%all_species,this%all_centres)
     call system_clock(tf)
     if (this%verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',0,' PS in',(tf-ts)/rate,' s'
    endif

    ! If we have asked for atomic predictions, we need to convert to atomic
    ! power spectra
    if (this%atomic) then
     ! Here, rather than using the fact that we know what these quantities are,
     ! we will instead use the size function -- so as to make this agnostic of
     ! what we know (a good direction for future code)
     this%tot_natoms = 0
     do i=1,size(this%PS,1)
      this%tot_natoms = this%tot_natoms + frames%natoms(i)
     enddo
     backup_nframes = frames%nframes
     backup_natmax  = frames%natmax
     frames%nframes = this%tot_natoms
     frames%natmax = 1
     if (allocated(this%atname_at)) deallocate(this%atname_at)
     if (allocated(this%natoms_at)) deallocate(this%natoms_at)
     allocate(this%PS_atomic(this%tot_natoms,1,size(this%PS,3), &
     &     size(this%PS,4)),this%natoms_at(this%tot_natoms),this%atname_at(this%tot_natoms))
     this%PS_atomic(:,:,:,:) = 0.d0
     this%natoms_at(:) = 0
     this%atname_at(:) = ''
     ! Populate atomic power spectrum array from original power spectrum
     k = 1
     do i=1,size(this%PS,1)
      this%PS_atomic(k:k+frames%natoms(i),1,:,:) = this%PS(i,1:frames%natoms(i),:,:)
      this%natoms_at(k:k+frames%natoms(i)-1) = frames%natoms(i)
      this%atname_at(k:k+frames%natoms(i)-1) = frames%atname(i,1:frames%natoms(i))
      k = k + frames%natoms(i)
     enddo
     ! Create new power spectrum array with the corrected shape
     deallocate(this%PS)
     allocate(this%PS(this%tot_natoms,1,size(this%PS_atomic,3),size(this%PS_atomic,4)))
     this%PS(:,:,:,:) = this%PS_atomic(:,:,:,:)
     deallocate(this%PS_atomic)
     if (this%do_scalar) then
      ! Repeat for scalar power spectrum if appropriate
      allocate(this%PS0_atomic(this%tot_natoms,1,size(this%PS0,3),size(this%PS0,4)))
      k = 1
      do i=1,size(this%PS0,1)
       this%PS0_atomic(k:k+frames%natoms(i),1,:,:) = this%PS0(i,1:frames%natoms(i),:,:)
       k = k + frames%natoms(i)
      enddo
      deallocate(this%PS0)
      allocate(this%PS0(this%tot_natoms,1,size(this%PS0_atomic,3),size(this%PS0_atomic,4)))
      this%PS0(:,:,:,:) = this%PS0_atomic(:,:,:,:)
      deallocate(this%PS0_atomic)
     endif
     ! Store number-of-atoms array to put it back later
     allocate(natoms_backup(size(frames%natoms)))
     natoms_backup = frames%natoms
     ! Create new number-of-atoms array with the corrected shape
     deallocate(frames%natoms)
     allocate(frames%natoms(this%tot_natoms))
     frames%natoms(:) = 1
    endif

    ! Get kernel
    if (allocated(this%natoms_tr)) deallocate(this%natoms_tr)
    allocate(this%natoms_tr(this%nmol))
    this%natoms_tr(:) = 1.d0
    call system_clock(ts)
    if (.not.this%do_scalar) then
     if (this%lm.eq.0) then
      this%ker = do_scalar_kernel(real(this%PS),this%PS_tr_lam,frames%nframes,this%nmol, &
     &     this%nfeat,this%nfeat,frames%natmax,1,dfloat(frames%natoms),this%natoms_tr,this%zeta,.false.)
     else
      if (this%zeta.gt.1) stop 'ERROR: zeta>1 has been selected but no scalar power spectrum given!'
      this%ker_lm = do_linear_spherical_kernel(real(this%PS),this%PS_tr_lam,frames%nframes,this%nmol, &
       &     this%nfeat,this%nfeat, &
       &     frames%natmax,1,dfloat(frames%natoms),this%natoms_tr,this%zeta,.false.,this%degen)
     endif
    else
     if (this%lm.eq.0) stop 'ERROR: you have asked for a spherical power spectrum but provided lambda=0!'
     this%ker_lm = do_nonlinear_spherical_kernel(real(this%PS),this%PS_tr_lam,real(this%PS0), &
     &     this%PS_tr_0,frames%nframes,this%nmol,this%nfeat,this%nfeat, &
     &     this%nfeat0,this%nfeat0,frames%natmax,1,dfloat(frames%natoms),this%natoms_tr,this%zeta,.false.,this%degen)
    endif
    call system_clock(tf)
    if (this%verbose) write(*,'(A,F6.3,A)') 'Got kernel in ',(tf-ts)/rate,' s'

    ! Get predictions
    if (allocated(this%prediction_lm_c)) deallocate(this%prediction_lm_c)
    allocate(this%prediction_lm_c(this%nmol,this%degen,this%nw))
    if (this%lm.eq.0) then
     this%prediction_lm_c = do_prediction_c(this%ker,this%wt_c,this%meanval_c,this%degen,frames%nframes, &
     &     this%nmol,this%nw)
    else
     this%prediction_lm_c = do_prediction_c(this%ker_lm,this%wt_c,this%meanval_c,this%degen,frames%nframes, &
     &     this%nmol,this%nw)
    endif

    if (this%atomic) then
     ! Rescale predictions by the number of atoms, so that they are related to
     ! the frame predictions by a summation
     do i=1,size(this%prediction_lm_c,1)
      this%prediction_lm_c(i,:,:) = this%prediction_lm_c(i,:,:) / this%natoms_at(i)
     enddo
     ! Set number of atoms equal to backed up version, in case we are applying
     ! multiple models
     deallocate(frames%natoms)
     allocate(frames%natoms(size(natoms_backup)))
     frames%natoms = natoms_backup
     deallocate(natoms_backup)
     frames%nframes = backup_nframes
     frames%natmax  = backup_natmax
    endif

end subroutine

!****************************************************************************************************************

subroutine rescale_predictions(this)
implicit none

    type(SAGPR_Model), intent(inout) :: this
    integer i,j

    ! Rescale committee predictions about their mean
    allocate(this%prediction_lm(size(this%prediction_lm_c,1),this%degen))
    this%prediction_lm(:,:) = 0.d0
    do i=1,this%degen
     do j=1,this%nw
      this%prediction_lm(:,i) = this%prediction_lm(:,i) + this%prediction_lm_c(:,i,j)
     enddo
     this%prediction_lm(:,i) = this%prediction_lm(:,i) / float(this%nw)
    enddo
    do j=1,this%nw
     this%prediction_lm_c(:,:,j) = this%prediction_lm(:,:) + (this%nu * (this%prediction_lm_c(:,:,j)-this%prediction_lm(:,:)))
    enddo
    deallocate(this%prediction_lm)

end subroutine

!****************************************************************************************************************

end module
