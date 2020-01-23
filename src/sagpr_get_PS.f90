program sagpr_get_PS
    use sagpr
    implicit none

    integer nargs,numkeys,ios,nlines,nframes,frnum,natmax,reals,bytes
    character(len=100), allocatable :: arg(:),keylist(:)
    integer lm,i,j,k,ii
    integer nmax,lmax,ncut
    real*8 rcut,sg,rs(3)
    character(len=100) ofile,sparse,fname
    logical periodic,readnext
    logical all_species(nelements),all_centres(nelements)
    real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
    character(len=4), allocatable :: atname(:,:)
    integer, allocatable :: natoms(:)
    character(len=1000), allocatable :: comment(:)
    character(len=1000) c1
    complex*16, allocatable, target :: sparse_data(:)
    complex*16, allocatable :: sparsification(:,:,:)
    character(len=3) symbol

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    lm = 0
    nmax = 8
    lmax = 6
    rcut = 4.d0
    sg = 0.3d0
    all_centres(:) = .false.
    all_species(:) = .false.
    sparse = ''
    fname = ''
    rs = (/0.d0,0.d0,0.d0/)
    ofile = ''
    periodic = .false.
    numkeys = 13
    allocate(keylist(numkeys))
    keylist = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-sf ','-f  ','-rs ','-o  ','-p  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
     if (arg(i).eq.'-n') read(arg(i+1),*) nmax
     if (arg(i).eq.'-l') read(arg(i+1),*) lmax
     if (arg(i).eq.'-rc') read(arg(i+1),*) rcut
     if (arg(i).eq.'-sg') read(arg(i+1),*) sg
     if (arg(i).eq.'-c') then
      readnext = .true.
      do k=i+1,nargs
       do j=1,numkeys
        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keylist(j))))
       enddo
       if (readnext) then
        read(arg(k),*) symbol
        do ii=1,nelements
         if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_centres(ii) = .true.
        enddo
       endif
      enddo
     endif
     if (arg(i).eq.'-s') then
      readnext = .true.
      do k=i+1,nargs
       do j=1,numkeys
        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keylist(j))))
       enddo
       if (readnext) then
        read(arg(k),*) symbol
        do ii=1,nelements
         if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_species(ii) = .true.
        enddo
       endif
      enddo
     endif
     if (arg(i).eq.'-sf') read(arg(i+1),*) sparse
     if (arg(i).eq.'-f') read(arg(i+1),*) fname
     if (arg(i).eq.'-rs') then
      read(arg(i+1),*) rs(1)
      read(arg(i+2),*) rs(2)
      read(arg(i+3),*) rs(3)
     endif
     if (arg(i).eq.'-o') read(arg(i+1),*) ofile
     if (arg(i).eq.'-p') periodic=.true.
    enddo

    ! Check for arguments that are required
    if (fname.eq.'') stop 'ERROR: filename required!'
    if (ofile.eq.'') stop 'ERROR: output file required!'
    if (sparse.eq.'') stop 'ERROR: sparsification file required!'

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

    ! If necessary, get sparsification details
    if (sparse.eq.'') then
     allocate(sparsification(1,1,1))
     sparsification(1,1,1) = -1.d0
     ncut = -1
    else
     open(unit=32,file=sparse,status='old',access='stream',form='unformatted')
     inquire(unit=32,size=bytes)
     reals = bytes / 16
     allocate(sparse_data(reals))
     read(32,pos=1) sparse_data
     close(32)
     ncut = sparse_data(1)
     if (reals.ne.(1 + (ncut*(ncut+1)))) then
      write(*,*) 'With ',reals,'elements instead of ',1 + (ncut*(ncut+1)),':'
      stop 'ERROR: incorrect number of elements in sparse file!'
     endif
     allocate(sparsification(2,ncut,ncut))
     do i=1,ncut
      sparsification(1,i,1) = sparse_data(i+1)
     enddo
     do i=1,ncut
      do j=1,ncut
       sparsification(2,i,j) = sparse_data(1 + ncut + (i-1)*ncut + j)
      enddo
     enddo
     deallocate(sparse_data)
    endif

    ! Get power spectrum
    call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut,sparsification,rs,periodic)

    ! Print power spectrum
    open(unit=33,file=ofile,access='stream',form='unformatted')
    write(33,pos=1) real(PS)
    close(33)

    ! Array deallocation
    deallocate(xyz,atname,natoms,comment,sparsification,cell,PS)

end program
