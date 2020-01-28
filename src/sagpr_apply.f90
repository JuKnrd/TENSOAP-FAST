program sagpr_apply
    use sagpr
    implicit none

    integer nargs,numkeys,ios,nlines,nframes,frnum,natmax,reals,bytes
    character(len=100), allocatable :: arg(:),keylist(:)
    integer lm,i,j,k,l,ii,nfeat,degen,nmol,nenv
    integer nmax,lmax,ncut,zeta,nfeat0,ncut0
    real*8 rcut,sg,rs(3),meanval
    character(len=100) ofile,sparse,fname,weights,ptrain
    logical periodic,readnext
    logical all_species(nelements),all_centres(nelements)
    real*8, allocatable :: xyz(:,:,:),cell(:,:,:)
    character(len=4), allocatable :: atname(:,:)
    integer, allocatable :: natoms(:)
    character(len=1000), allocatable :: comment(:)
    character(len=1000) c1
    real*8, allocatable, target :: raw_PS(:),all_wt(:)
    real*8, allocatable :: PS_tr_lam(:,:,:),PS_tr_0(:,:,:),ker(:,:), ker_lm(:,:,:,:),natoms_tr(:)
    real*8, allocatable :: prediction_lm(:,:),wt(:)
    complex*16, allocatable, target :: sparse_data(:)
    complex*16, allocatable :: sparsification(:,:,:),sparsification0(:,:,:)
    character(len=3) symbol
    logical do_scalar

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
    weights = ''
    rs = (/0.d0,0.d0,0.d0/)
    ofile = 'prediction.out'
    periodic = .false.
    zeta = 1
    ptrain = ''
    numkeys = 16
    allocate(keylist(numkeys))
    keylist = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-sf ','-f  ','-rs ','-o  ','-p  ','-z  ','-w  ','-pt ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
     if (arg(i).eq.'-n') read(arg(i+1),*) nmax
     if (arg(i).eq.'-l') read(arg(i+1),*) lmax
     if (arg(i).eq.'-rc') read(arg(i+1),*) rcut
     if (arg(i).eq.'-sg') read(arg(i+1),*) sg
     if (arg(i).eq.'-z') read(arg(i+1),*) zeta
     if (arg(i).eq.'-w') read(arg(i+1),*) weights
     if (arg(i).eq.'-pt') read(arg(i+1),*) ptrain
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
    if (weights.eq.'') stop 'ERROR: weights file required!'
    if (sparse.eq.'') stop 'ERROR: sparsification file required!'
    if (ptrain.eq.'') stop 'ERROR: training power spectra required!'

    ! Read in power spectrum file(s)
    degen = 2*lm + 1
    open(unit=41,file=ptrain,status='old',access='stream',form='unformatted')
    inquire(unit=41,size=bytes)
    reals = bytes/8
    allocate(raw_PS(reals))
    read(41,pos=1) raw_PS
    do_scalar = (raw_PS(1).ne.0.d0)
    nmol = int(raw_PS(2))
    nfeat = int(raw_PS(3))
    i = 3
    if (do_scalar) then
     nfeat0 = int(raw_PS(4))
     i = 4
    endif
    allocate(PS_tr_lam(nmol,degen,nfeat))
    do j=1,nmol
     do k=1,degen
      do l=1,nfeat
       i = i + 1
       PS_tr_lam(j,k,l) = raw_PS(i)
      enddo
     enddo
    enddo
    if (do_scalar) then
     allocate(PS_tr_0(nmol,1,nfeat))
     do j=1,nmol
      do k=1,nfeat0
       PS_tr_0(j,1,k) = raw_PS(i)
      enddo
     enddo
    endif
    if (i.ne.reals) stop 'ERROR: different file size to that expected for power spectrum!'

    ! Get sparsification details
    open(unit=32,file=sparse,status='old',access='stream',form='unformatted')
    inquire(unit=32,size=bytes)
    reals = bytes / 16
    allocate(sparse_data(reals))
    read(32,pos=1) sparse_data
    close(32)
    ncut = int(sparse_data(1))
    allocate(sparsification(2,ncut,ncut))
    i = 1
    do j=1,ncut
     i = i + 1
     sparsification(1,j,1) = sparse_data(i)
    enddo
    do j=1,ncut
     do k=1,ncut
      i = i + 1
      sparsification(2,j,k) = sparse_data(i)
     enddo
    enddo
    if (do_scalar) then
     i = i + 1
     ncut0 = int(sparse_data(i))
     allocate(sparsification0(2,ncut0,ncut0))
     do j=1,ncut0
      i = i + 1
      sparsification0(1,j,1) = sparse_data(i)
     enddo
     do j=1,ncut0
      do k=1,ncut0
       i = i + 1
       sparsification0(2,j,k) = sparse_data(i)
      enddo
     enddo
    endif
    if (i.ne.reals) stop 'ERROR: different file size to that expected for sparsification!'
    deallocate(sparse_data)

    if (ncut.ne.nfeat) stop 'ERROR: ncut .ne. nfeat!'
    if (do_scalar) then
     if (ncut0.ne.nfeat0) stop 'ERROR: ncut0 .ne. nfeat0!'
    endif

    ! Get weights
    open(unit=31,file=weights,status='old',access='stream',form='unformatted')
    inquire(unit=31,size=bytes)
    reals = bytes/8
    allocate(all_wt(reals),wt(reals-1))
    read(31,pos=1) all_wt
    meanval = all_wt(1)
    do i=2,reals
     wt(i-1) = all_wt(i)
    enddo
    nenv = (reals-1) / degen
    if (nenv .ne. nmol) stop 'ERROR: size of weight vector does not match number of environments!'
    close(31)

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

    ! Get power spectrum
    if (.not.do_scalar) then
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut,sparsification,rs,periodic)
    else
     stop 'NOT YET SET UP FOR DOING ZETA>1'
    endif

    ! Get kernel
    allocate(natoms_tr(nmol))
    natoms_tr(:) = 1.d0
    if (.not.do_scalar) then
     if (lm.eq.0) then
      ker = do_scalar_kernel(real(PS),PS_tr_lam,nframes,nmol,nfeat,nfeat,natmax,1,dfloat(natoms),natoms_tr,zeta,.false.)
     else
      if (zeta.gt.1) stop 'ERROR: zeta>1 has been selected but no scalar power spectrum given!'
      ker_lm = do_linear_spherical_kernel(real(PS),PS_tr_lam,nframes,nmol,nfeat,nfeat, &
     &     natmax,1,dfloat(natoms),natoms_tr,zeta,.false.,degen)
     endif
    else
     stop 'NOT YET SET UP FOR DOING ZETA>1'
    endif

    ! Get predictions
    allocate(prediction_lm(nmol,degen))
    prediction_lm = do_prediction(ker,wt,meanval,degen,nframes,nmol)

    ! Print predictions
    open(unit=33,file=ofile)
    do i=1,nframes
     write(33,*) (prediction_lm(i,j),j=1,degen)
    enddo
    close(33)

    ! Array deallocation
    deallocate(xyz,atname,natoms,comment,sparsification,cell,PS,PS_tr_lam,natoms_tr)
    if (allocated(PS_tr_0)) deallocate(PS_tr_0)
    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)

end program
