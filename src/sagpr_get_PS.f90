program sagpr_predict
    use sagpr
    implicit none

    integer nargs,numkeys
    character(len=100), allocatable :: arg(:),keylist(:)
    integer lm,i,j,nmol,reals,bytes,degen,nenv,mu,nu,k,ii
    integer nmax,lmax,ncut
    real*8 rcut,sg,rs(3)
    character(len=100) weights,ofile,kerfile,sparse,fname
    real*8, allocatable, target :: all_wt(:),raw_ker(:)
    real*8, allocatable :: wt(:),prediction_lm(:,:),ker(:,:,:,:)
    logical periodic,readnext
    real*8 meanval
    integer, parameter :: nelements = 200
    logical all_species(nelements),all_centres(nelements)

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

    ! Array deallocation
    !TODO: deallocation

end program
