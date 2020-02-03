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
 logical all_species(nelements),all_centres(nelements)

 contains

!****************************************************************************************************************
subroutine set_defaults()
 implicit none

  lm = 0
  nmax = 8
  lmax = 6
  rcut = 4.d0
  sg = 0.3d0
    all_centres(:) = .false.
    all_species(:) = .false.
  rs = (/0.d0,0.d0,0.d0/)
  periodic = .false.
  zeta = 1

end subroutine

!****************************************************************************************************************
function process_hyperparameters(model)
 implicit none

 character(len=100) model,args
 character(len=100), allocatable :: process_hyperparameters(:),keys2(:)
 integer i,j,k,ii,numkeys
 logical new_arg,readnext
 character(len=3) symbol

  ! Read in hyperparameters
  open(unit=21,file=trim(adjustl(model))//'.hyp',status='old')
  read(21,'(A)') args
  nargs = 0
  new_arg = .false.
  do i=1,len(trim(adjustl(args)))
   if (args(i:i).ne.' ') then
    if (.not.new_arg) then
     new_arg = .true.
     nargs = nargs + 1
    endif
   else
    ! We have reached the end of an argument
    new_arg = .false.
   endif
  enddo
  allocate(process_hyperparameters(nargs+1))
  read(args,*) (process_hyperparameters(i),i=1,nargs)
  process_hyperparameters(nargs+1) = 'NULL'

  ! Get hyperparameters from model file
  numkeys = 11
  allocate(keys2(numkeys))
  keys2 = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-rs ','-p  ','-z  ','NULL'/)
  do i=1,nargs
   process_hyperparameters(i) = trim(adjustl(process_hyperparameters(i)))
   if (process_hyperparameters(i).eq.'-lm') read(process_hyperparameters(i+1),*) lm
   if (process_hyperparameters(i).eq.'-n') read(process_hyperparameters(i+1),*) nmax
   if (process_hyperparameters(i).eq.'-l') read(process_hyperparameters(i+1),*) lmax
   if (process_hyperparameters(i).eq.'-rc') read(process_hyperparameters(i+1),*) rcut
   if (process_hyperparameters(i).eq.'-sg') read(process_hyperparameters(i+1),*) sg
   if (process_hyperparameters(i).eq.'-z') read(process_hyperparameters(i+1),*) zeta
   if (process_hyperparameters(i).eq.'-c') then
    readnext = .true.
    do k=i+1,nargs
     do j=1,numkeys
      readnext = readnext.and.(trim(adjustl(process_hyperparameters(k))).ne.trim(adjustl(keys2(j))))
     enddo
     if (readnext) then
      read(process_hyperparameters(k),*) symbol
      do ii=1,nelements
       if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_centres(ii) = .true.
      enddo
     endif
    enddo
   endif
   if (process_hyperparameters(i).eq.'-s') then
    readnext = .true.
    do k=i+1,nargs
     do j=1,numkeys
      readnext = readnext.and.(trim(adjustl(process_hyperparameters(k))).ne.trim(adjustl(keys2(j))))
     enddo
     if (readnext) then
      read(process_hyperparameters(k),*) symbol
      do ii=1,nelements
       if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_species(ii) = .true.
      enddo
     endif
    enddo
   endif
   if (process_hyperparameters(i).eq.'-rs') then
    read(process_hyperparameters(i+1),*) rs(1)
    read(process_hyperparameters(i+2),*) rs(2)
    read(process_hyperparameters(i+3),*) rs(3)
   endif
   if (process_hyperparameters(i).eq.'-p') periodic=.true.
  enddo
  ! Set degeneracy
  degen = 2*lm + 1

end function

!****************************************************************************************************************
subroutine get_model(model)
 implicit none

  character(len=100) model
  real*8, allocatable, target :: raw_model(:)
  real*8 a1,a2
  integer i,j,k,l

   ! Read in power spectrum file(s)
   open(unit=41,file=trim(adjustl(model))//'.mdl',status='old',access='stream',form='unformatted')
   inquire(unit=41,size=bytes)
   reals = bytes/8
   allocate(raw_model(reals))
   read(41,pos=1) raw_model
   do_scalar = (raw_model(1).ne.0.d0)
   nmol = int(raw_model(2))
   nfeat = int(raw_model(3))
   i = 3
   if (do_scalar) then
    nfeat0 = int(raw_model(4))
    i = 4
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
   i = i + 1
   meanval = raw_model(i)
   allocate(wt(nmol*degen))
   do j=1,nmol*degen
    i = i + 1
    wt(j) = raw_model(i)
   enddo
   if (i.ne.reals) stop 'ERROR: diferent file size to that expected for model!'

end subroutine

!****************************************************************************************************************

subroutine read_fifo(un,periodic)
 implicit none

  integer ios,un,nat
  character(len=1000) line,c1
  integer i,j,ii
  logical periodic

   nframes = 1
   read(un,'(A)') line
   read(line,*,iostat=ios) nat
   if (ios.ne.0) then
    call system('rm my_fifo_in my_fifo_out')
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
  write(*,*) 'Got input frame with ',nat,'atoms'

end subroutine

!****************************************************************************************************************

subroutine read_xyz(fname,periodic)
 implicit none

  character(len=100) fname
  character(len=1000) c1
  integer ios,i,j,frnum,k,ii
  logical periodic

    ! Read in XYZ file
    open(unit=31,file=fname,status='old')
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
      write(*,*) 'At frame ',k,':'
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

end module
