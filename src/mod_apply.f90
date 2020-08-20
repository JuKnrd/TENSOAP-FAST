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
 integer lm,nmax,lmax,zeta,degen
 real*8 rcut,sg,rs(3)
 logical periodic
 integer nmax0,lmax0
 real*8 rcut0,sg0,rs0(3)
 ! Defaults
 integer, parameter :: nmax_default = 8, lmax_default = 6
 real*8, parameter :: rcut_default = 4.d0,sg_default = 0.3d0,rs_default(3) = (/0.d0,0.d0,0.d0/)
 ! Committee models
 logical committee
 integer nw
 real*8, allocatable :: wt_c(:,:),meanval_c(:),prediction_lm_c(:,:,:)
 real*8 nu
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
subroutine get_model(model)
 implicit none

  character(len=100) model
  real*8, allocatable, target :: raw_model(:)
  real*8 a1,a2
  integer i,j,k,l,ncen,nspec

   ! Read in power spectrum file(s)
   open(unit=41,file=trim(adjustl(model)),status='old',access='stream',form='unformatted')
   inquire(unit=41,size=bytes)
   reals = bytes/8
   allocate(raw_model(reals))
   read(41,pos=1) raw_model
   do_scalar = (raw_model(1).ne.0.d0)
   lm = int(raw_model(2))
   degen=2*lm + 1
   zeta = int(raw_model(3))
   periodic = (int(raw_model(4)).eq.1)
   nmol = int(raw_model(5))
   nfeat = int(raw_model(6))
   i = 6
   if (do_scalar) then
    nfeat0 = int(raw_model(7))
    i = 7
   endif
   if (committee) then
    i = i + 1
    nw = int(raw_model(i))
   endif
   allocate(PS_tr_lam(nmol,1,degen,nfeat))
   do j=1,nmol
    do k=1,degen
     do l=1,nfeat
      i = i + 1
      PS_tr_lam(j,1,k,l) = raw_model(i)
     enddo
    enddo
   enddo
   if (do_scalar) then
    allocate(PS_tr_0(nmol,1,1,nfeat))
    do j=1,nmol
     do k=1,nfeat0
      i = i + 1
      PS_tr_0(j,1,1,k) = raw_model(i)
     enddo
    enddo
   endif

   ! Get sparsification details
   i = i + 1
   ncut = int(raw_model(i))
   allocate(sparsification(2,ncut,ncut))
   do j=1,ncut
    i = i + 1
    sparsification(1,j,1) = raw_model(i)
   enddo
   do j=1,ncut
    do k=1,ncut
     i = i + 1
     a1 = raw_model(i)
     i = i + 1
     a2 = raw_model(i)
     sparsification(2,j,k) = dcmplx(a1,a2)
    enddo
   enddo
   if (do_scalar) then
    i = i + 1
    ncut0 = int(raw_model(i))
    allocate(sparsification0(2,ncut0,ncut0))
    do j=1,ncut0
     i = i + 1
     sparsification0(1,j,1) = raw_model(i)
    enddo
    do j=1,ncut0
     do k=1,ncut0
      i = i + 1
      a1 = raw_model(i)
      i = i + 1
      a2 = raw_model(i)
      sparsification0(2,j,k) = dcmplx(a1,a2)
     enddo
    enddo
   endif

   if (ncut.ne.nfeat) stop 'ERROR: ncut .ne. nfeat!'
   if (do_scalar) then
    if (ncut0.ne.nfeat0) stop 'ERROR: ncut0 .ne. nfeat0!'
   endif

   ! Get weights
   if (committee) then
    allocate(meanval_c(nw),wt_c(nmol*degen,nw))
    do k=1,nw
     i = i + 1
     meanval_c(k) = raw_model(i)
    enddo
    do k=1,nw
     do j=1,nmol*degen
      i = i + 1
      wt_c(j,k) = raw_model(i)
     enddo
    enddo
   else
    i = i + 1
    meanval = raw_model(i)
    allocate(wt(nmol*degen))
    do j=1,nmol*degen
     i = i + 1
     wt(j) = raw_model(i)
    enddo
   endif

   ! Get hyperparameters
   i = i + 1
   nmax = int(raw_model(i))
   if (nmax.eq.-1) nmax=nmax_default
   i = i + 1
   lmax = int(raw_model(i))
   if (lmax.eq.-1) lmax=lmax_default
   i = i + 1
   rcut = raw_model(i)
   if (rcut.lt.0.d0) rcut=rcut_default
   i = i + 1
   sg = raw_model(i)
   if (sg.lt.0.d0) sg=sg_default
   i = i + 1
   rs(1) = raw_model(i)
   i = i + 1
   rs(2) = raw_model(i)
   i = i + 1
   rs(3) = raw_model(i)
   if (rs(3).eq.0.d0) rs=rs_default
   if (do_scalar) then
    i = i + 1
    nmax0 = int(raw_model(i))
    if (nmax0.eq.-1) nmax0=nmax_default
    i = i + 1
    lmax0 = int(raw_model(i))
    if (lmax0.eq.-1) lmax0=lmax_default
    i = i + 1
    rcut0 = raw_model(i)
    if (rcut0.lt.0.d0) rcut0=rcut_default
    i = i + 1
    sg0 = raw_model(i)
    if (sg0.lt.0.d0) sg0=sg_default
    i = i + 1
    rs0(1) = raw_model(i)
    i = i + 1
    rs0(2) = raw_model(i)
    i = i + 1
    rs0(3) = raw_model(i)
    if (rs0(3).eq.0.d0) rs0=rs_default
   endif
   i = i + 1
   ncen = int(raw_model(i))
   if (ncen.gt.0) then
    do j=1,ncen
     i = i + 1
     all_centres(int(raw_model(i))) = .true.
    enddo
   endif
   i = i + 1
   nspec = int(raw_model(i))
   if (nspec.gt.0) then
    do j=1,ncen
     i = i + 1
     all_species(int(raw_model(i))) = .true.
    enddo
   endif
   if (committee) then
    i = i + 1
    nu = raw_model(i)**0.5
   endif

   if (i.ne.reals) stop 'ERROR: different file size to that expected for model!'

end subroutine

!****************************************************************************************************************

subroutine read_fifo(un,periodic,model)
 implicit none

  integer ios,un,nat
  character(len=1000) line,c1
  character(len=100) model
  integer i,j,ii
  logical periodic

   nframes = 1
   read(un,'(A)') line
   read(line,*,iostat=ios) nat
   if (ios.ne.0) then
    call system('rm '//trim(adjustl(model))//'.{in,out}')
    write(*,*) 'Non-standard input detected; now stopping'
    stop
   endif
   natmax = nat
   if (allocated(xyz)) deallocate(xyz,atname,natoms,comment)
   allocate(xyz(nframes,natmax,3),atname(nframes,natmax),natoms(nframes),comment(nframes))
   natoms(1) = nat
   read(un,'(A)') comment(1)
   do j=1,natoms(1)
    read(un,*) atname(1,j),(xyz(1,j,ii),ii=1,3)
   enddo

  ! Get cell data
  if (allocated(cell)) deallocate(cell)
  allocate(cell(nframes,3,3))
  if (.not.periodic) then
   cell(:,:,:) = 0.d0
  else
   do i=1,nframes
    ios = index(comment(i),'Lattice')
    if (ios.eq.0) stop 'ERROR: input file is not periodic!'
    c1 = comment(i)
    c1 = c1(ios:len(c1))
    ios = index(c1,'"')
    c1 = c1(ios+1:len(c1))
    ios = index(c1,'"')
    c1 = c1(1:ios-1)
    read(c1,*) cell(i,1,1),cell(i,1,2),cell(i,1,3),cell(i,2,1),cell(i,2,2),cell(i,2,3),cell(i,3,1),cell(i,3,2),cell(i,3,3)
   enddo
  endif
  if (verbose) write(*,*) 'Got input frame with ',nat,'atoms'

end subroutine

!****************************************************************************************************************

subroutine read_frame(un,periodic)
 implicit none

  integer ios,un,nat
  character(len=1000) line,c1
  character(len=100) model
  integer i,j,ii
  logical periodic

   nframes = 1
   read(un,'(A)',iostat=ios) line
   if (ios.ne.0) then
    if (verbose) write(*,*) 'End-of-file detected'
    stop
   endif
   read(line,*,iostat=ios) nat
   if (ios.ne.0) then
    call system('rm '//trim(adjustl(model))//'.{in,out}')
    write(*,*) 'Non-standard input detected; now stopping'
    stop
   endif
   natmax = nat
   if (allocated(xyz)) deallocate(xyz,atname,natoms,comment)
   allocate(xyz(nframes,natmax,3),atname(nframes,natmax),natoms(nframes),comment(nframes))
   natoms(1) = nat
   read(un,'(A)') comment(1)
   do j=1,natoms(1)
    read(un,*) atname(1,j),(xyz(1,j,ii),ii=1,3)
   enddo

  ! Get cell data
  if (allocated(cell)) deallocate(cell)
  allocate(cell(nframes,3,3))
  if (.not.periodic) then
   cell(:,:,:) = 0.d0
  else
   do i=1,nframes
    ios = index(comment(i),'Lattice')
    if (ios.eq.0) stop 'ERROR: input file is not periodic!'
    c1 = comment(i)
    c1 = c1(ios:len(c1))
    ios = index(c1,'"')
    c1 = c1(ios+1:len(c1))
    ios = index(c1,'"')
    c1 = c1(1:ios-1)
    read(c1,*) cell(i,1,1),cell(i,1,2),cell(i,1,3),cell(i,2,1),cell(i,2,2),cell(i,2,3),cell(i,3,1),cell(i,3,2),cell(i,3,3)
   enddo
  endif
  if (verbose) write(*,*) 'Got input frame with ',nat,'atoms'

end subroutine

!****************************************************************************************************************

subroutine read_xyz(fname,periodic)
 implicit none

  character(len=100) fname
  character(len=1000) c1
  integer ios,i,j,frnum,k,ii
  logical periodic

    ! Read in XYZ file
    open(unit=31,file=fname,status='old',iostat=ios)
    if (ios.ne.0) stop 'ERROR: input file does not exist!'
    nlines = 0
    do
     read(31,*,iostat=ios)
     if (ios.ne.0) exit
     nlines = nlines + 1
    enddo
    close(31)

    open(unit=31,file=fname,status='old')
    nframes = 0
    natmax = 0
    do i=1,nlines
     read(31,*,iostat=ios) frnum
     if (ios.eq.0) then
      nframes = nframes + 1
      natmax = max(natmax,frnum)
     endif
    enddo
    close(31)

    if (allocated(xyz)) deallocate(xyz,atname,natoms,comment)
    allocate(xyz(nframes,natmax,3),atname(nframes,natmax),natoms(nframes),comment(nframes))
    xyz(:,:,:) = 0.d0
    atname(:,:) = ''
    open(unit=31,file=fname,status='old')
    k = 0
    do i=1,nframes
     read(31,*,iostat=ios) frnum
     if (ios.eq.0) then
      k = k + 1
      natoms(k) = frnum
      read(31,'(A)') comment(k)
      do j=1,natoms(k)
       read(31,*) atname(k,j),(xyz(k,j,ii),ii=1,3)
      enddo
     else
      if (verbose) write(*,*) 'At frame ',k,':'
      stop 'ERROR: there should be an atom number here!'
     endif
    enddo
    close(31)

    ! Get cell data
    if (allocated(cell)) deallocate(cell)
    allocate(cell(nframes,3,3))
    if (.not.periodic) then
     cell(:,:,:) = 0.d0
    else
     do i=1,nframes
      ios = index(comment(i),'Lattice')
      if (ios.eq.0) stop 'ERROR: input file is not periodic!'
      c1 = comment(i)
      c1 = c1(ios:len(c1))
      ios = index(c1,'"')
      c1 = c1(ios+1:len(c1))
      ios = index(c1,'"')
      c1 = c1(1:ios-1)
      read(c1,*) cell(i,1,1),cell(i,1,2),cell(i,1,3),cell(i,2,1),cell(i,2,2),cell(i,2,3),cell(i,3,1),cell(i,3,2),cell(i,3,3)
     enddo
    endif

    return

end subroutine
!****************************************************************************************************************
subroutine predict_frame(rate)
 implicit none

    integer ts,tf
    real*8 rate

    ! Get power spectrum
    if (.not.do_scalar) then
     call system_clock(ts)
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,F6.3,A)') 'Got PS in',(tf-ts)/rate,' s'
    else
     call system_clock(ts)
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',lm,' PS in',(tf-ts)/rate,' s'
     call system_clock(ts)
     call do_power_spectrum_scalar(xyz,atname,natoms,cell,nframes,natmax,0,nmax0,lmax0,rcut0,sg0, &
     &     ncut0,sparsification0,rs0,periodic,.true.)
     call system_clock(tf)
     if (verbose) write(*,'(A,I2,A,F6.3,A)') 'Got L=',0,' PS in',(tf-ts)/rate,' s'
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
    if (committee) then
     if (allocated(prediction_lm_c)) deallocate(prediction_lm_c)
     allocate(prediction_lm_c(nmol,degen,nw))
     if (lm.eq.0) then
      prediction_lm_c = do_prediction_c(ker,wt_c,meanval_c,degen,nframes,nmol,nw)
     else
      prediction_lm_c = do_prediction_c(ker_lm,wt_c,meanval_c,degen,nframes,nmol,nw)
     endif
    else
     if (allocated(prediction_lm)) deallocate(prediction_lm)
     allocate(prediction_lm(nmol,degen))
     if (lm.eq.0) then
      prediction_lm = do_prediction(ker,wt,meanval,degen,nframes,nmol)
     else
      prediction_lm = do_prediction(ker_lm,wt,meanval,degen,nframes,nmol)
     endif
    endif

end subroutine

!****************************************************************************************************************

end module
