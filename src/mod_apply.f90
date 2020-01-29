module apply
 use sagpr

 ! Variables for getting coordinates
 real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
 character(len=4), allocatable :: atname(:,:)
 character(len=1000), allocatable :: comment(:)
 integer, allocatable :: natoms(:)
 integer nframes,nlines,natmax
 ! Variables for getting models
 integer nargs
 real*8, allocatable :: PS_tr_lam(:,:,:,:),PS_tr_0(:,:,:,:),ker(:,:), ker_lm(:,:,:,:),natoms_tr(:)
 real*8, allocatable :: prediction_lm(:,:),wt(:)
 complex*16, allocatable :: sparsification(:,:,:),sparsification0(:,:,:)
 ! Parameters for PS and kernel building
 integer lm,nmax,lmax,zeta
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
 character(len=100), allocatable :: process_hyperparameters(:)
 integer i
 logical new_arg

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

end function

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

    ! Get cell data
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
