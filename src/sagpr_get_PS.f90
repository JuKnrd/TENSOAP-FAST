program sagpr_predict
    use sagpr
    implicit none

    integer nargs,numkeys,ios,nlines,nframes,frnum,natmax
    character(len=100), allocatable :: arg(:),keylist(:)
    integer lm,i,j,k,ii
    integer nmax,lmax,ncut
    real*8 rcut,sg,rs(3)
    character(len=100) ofile,sparse,fname
    logical periodic,readnext
    integer, parameter :: nelements = 200
    logical all_species(nelements),all_centres(nelements)
    real*8, allocatable :: xyz(:,:,:),PS(:,:,:,:)
    character(len=4), allocatable :: atname(:,:)
    integer, allocatable :: natoms(:)
    character(len=100), allocatable :: comment(:)

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
    ncut = -1
    sparse = ''
    fname = ''
    rs = (/0.d0,0.d0,0.d0/)
    ofile = ''
    periodic = .false.
    numkeys = 14
    allocate(keylist(numkeys))
    keylist = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-nc ','-sf ','-f  ','-rs ','-o  ','-p  ','NULL'/)
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
        read(arg(k),*) ii
        all_centres(ii) = .true.
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
        read(arg(k),*) ii
        all_species(ii) = .true.
       endif
      enddo
     endif
     if (arg(i).eq.'-nc') read(arg(i+1),*) ncut
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
      read(31,*) comment(k)
      do j=1,natoms(k)
       read(31,*) atname(k,j),(xyz(k,j,ii),ii=1,3)
      enddo
     else
      write(*,*) 'At frame ',k,':'
      stop 'ERROR: there should be an atom number here!'
     endif
    enddo

    ! Get power spectrum
    PS = do_power_spectrum()

    ! Print power spectrum
    open(unit=32,file=ofile,access='stream',form='unformatted')
    write(32,pos=1) PS
    close(32)

    ! Array deallocation
    deallocate(xyz,atname,natoms,comment,PS)

end program
