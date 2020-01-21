program sagpr_get_kernel
    use sagpr
    use iso_c_binding
    implicit none

    integer nargs,numkeys,ios
    character(len=100), allocatable :: arg(:),keylist(:)
    character(len=100) psname(2),psname0(2),ofile,scaling(2)
    integer zeta,lm,degen
    integer nunits,nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,nfeat0_1,nfeat0_2
    real*8, allocatable, target :: natoms_1(:),natoms_2(:), raw_PS_1(:),raw_PS_2(:),raw_PS0_1(:),raw_PS0_2(:)
    real*8, pointer :: PS_1(:,:,:),PS_2(:,:,:),PS0_1(:,:,:),PS0_2(:,:,:)
    real*8, pointer :: PS_1_lm(:,:,:,:),PS_2_lm(:,:,:,:)
    real*8, allocatable :: PS_1_fast(:,:),PS_2_fast(:,:),PS_1_fast_lm(:,:,:),PS_2_fast_lm(:,:,:)

    integer :: bytes,reals,natmax,nmol,nfeat,i,j,ii,jj,j0,mu,nu

    real*8, allocatable :: ker(:,:), ker_lm(:,:,:,:)
    real*8 k0
    logical :: hermiticity = .false.,readnext

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    zeta    = 1
    lm      = 0
    psname  = (/'',''/)
    psname0 = (/'',''/)
    scaling = (/'',''/)
    ofile   = ''
    numkeys = 7
    allocate(keylist(numkeys))
    keylist=(/'-z  ','-ps ','-ps0','-s  ','-o  ','-lm ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-z') read(arg(i+1),*) zeta
     if (arg(i).eq.'-ps') then
      read(arg(i+1),*) psname(1)
      readnext = .true.
      do j=1,numkeys
       readnext = readnext.and.(arg(i+2).ne.trim(adjustl(keylist(j))))
      enddo
      if (readnext) then
       read(arg(i+2),*) psname(2)
      else
       psname(2) = psname(1)
      endif
      if (psname(1).eq.psname(2)) hermiticity=.true.
     endif
     if (arg(i).eq.'-ps0') then
      read(arg(i+1),*) psname0(1)
      readnext = .true.
      do j=1,numkeys
       readnext = readnext.and.(arg(i+2).ne.trim(adjustl(keylist(j))))
      enddo
      if (readnext) then
       read(arg(i+2),*) psname0(2)
      else
       psname0(2) = psname0(1)
      endif
     endif
     if (arg(i).eq.'-s') then
      read(arg(i+1),*) scaling(1)
      readnext = .true.
      do j=1,numkeys
       readnext = readnext.and.(arg(i+2).ne.trim(adjustl(keylist(j))))
      enddo
      if (readnext) then
       read(arg(i+2),*) scaling(2)
      else
       scaling(2) = scaling(1)
      endif
     endif
     if (arg(i).eq.'-o') read(arg(i+1),*) ofile
     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
    enddo
    deallocate(keylist)

    ! Check for arguments that are required
    if (ofile.eq.'') stop 'ERROR: output file required!'
    if (psname(1).eq.'' .or. psname(2).eq.'') stop 'ERROR: power spectrum file(s) required!'
    if (scaling(1).eq.'' .or. scaling(2).eq.'') stop 'ERROR: scaling file(s) required!'

    ! Read in scaling file(s); if a file doesn't exist, then it is the number of 1s that has been entered
    do i=1,2
     open(unit=40,file=scaling(i),status='old',access='stream',form='unformatted',iostat=ios)
     if (ios.eq.0) then
      if (i.eq.1) then
       inquire(unit=40,size=bytes)
       reals = bytes/8
       allocate(natoms_1(reals))
       nmol_1 = reals
       read(40,pos=1) natoms_1
       close(40)
      else
       inquire(unit=40,size=bytes)
       reals = bytes/8
       allocate(natoms_2(reals))
       nmol_2 = reals
       read(40,pos=1) natoms_2
       close(40)
      endif
     else
      read(scaling(i),*) nunits
      if (i.eq.1) then
       allocate(natoms_1(nunits))
       natoms_1(:) = 1.d0
      else
       allocate(natoms_2(nunits))
       natoms_2(:) = 1.d0
      endif
     endif
    enddo

    ! Read in power spectrum file(s)
    degen = 2*lm + 1
    do i=1,2
     open(unit=41,file=psname(i),status='old',access='stream',form='unformatted')
     inquire(unit=41,size=bytes)
     reals = bytes/8
     if (i.eq.1) then
      allocate(raw_PS_1(reals))
      read(41,pos=1) raw_PS_1
      close(41)
      natmax_1 = maxval(natoms_1)
      nfeat_1 = reals / (nmol_1*natmax_1*degen)
      call C_F_POINTER(C_LOC(raw_PS_1),PS_1_lm,[nfeat_1,degen,natmax_1,nmol_1])
     else
      allocate(raw_PS_2(reals))
      read(41,pos=1) raw_PS_2
      close(41)
      natmax_2 = maxval(natoms_2)
      nfeat_2 = reals / (nmol_2*natmax_2*degen)
      call C_F_POINTER(C_LOC(raw_PS_2),PS_2_lm,[nfeat_2,degen,natmax_2,nmol_2])
     endif
    enddo

    ! Build kernels
    if (lm.eq.0) then

     ker = do_scalar_kernel(PS_1_lm,PS_2_lm,nmol_1,nmol_2,nfeat_1,nfeat_2,natmax_1,natmax_2,natoms_1,natoms_2,zeta,hermiticity)

!        ! Get kernels
!        if (zeta.eq.1) then
!         allocate(PS_1_fast(nfeat_1,nmol_1),PS_2_fast(nfeat_2,nmol_2))
!         PS_1_fast(:,:) = 0.d0
!         PS_2_fast(:,:) = 0.d0
!         do i=1,nmol_1
!          do ii=1,int(natoms_1(i))
!           PS_1_fast(:,i) = PS_1_fast(:,i) + PS_1_lm(:,1,ii,i) / natoms_1(i)
!          enddo
!         enddo
!         do i=1,nmol_2
!          do ii=1,int(natoms_2(i))
!           PS_2_fast(:,i) = PS_2_fast(:,i) + PS_2_lm(:,1,ii,i) / natoms_2(i)
!          enddo
!         enddo
!         allocate(ker(nmol_1,nmol_2))
!         ker(:,:) = 0.d0
!         do i=1,nmol_1
!          j0=1
!          if (hermiticity) j0=i
!          !$OMP PARALLEL DO SHARED(ker,PS_1_fast,PS_2_fast,natoms_1,natoms_2)
!          do j=j0,nmol_2
!           ker(i,j) = (dot_product(PS_1_fast(:,i),PS_2_fast(:,j)))
!           if (hermiticity) ker(j,i) = ker(i,j)
!          enddo
!          !$OMP END PARALLEL DO
!         enddo          
!        else
!         allocate(ker(nmol_1,nmol_2))
!         ker(:,:) = 0.d0
!         do i=1,nmol_1
!          j0=1
!          if (hermiticity) j0=i
!          !$OMP PARALLEL DO SHARED(ker,PS_1,PS_2,natoms_1,natoms_2) PRIVATE(ii,jj)
!          do j=j0,nmol_2
!           do ii=1,int(natoms_1(i))
!            do jj=1,int(natoms_2(j))
!             ker(i,j) = ker(i,j) + (dot_product(PS_1_lm(:,1,ii,i),PS_2_lm(:,1,jj,j)))**zeta
!            enddo
!           enddo
!           ker(i,j) = ker(i,j) / (natoms_1(i)*natoms_2(j))
!           if (hermiticity) ker(j,i) = ker(i,j)
!          enddo
!          !$OMP END PARALLEL DO
!         enddo
!        endif

        deallocate(arg,natoms_1,natoms_2,raw_PS_1,raw_PS_2)
        if (allocated(PS_1_fast)) deallocate(PS_1_fast)
        if (allocated(PS_2_fast)) deallocate(PS_2_fast)

    else
!
!        if (zeta.eq.1) then
!         ! Build kernel
!         allocate(ker_lm(nmol_1,nmol_2,degen,degen))
!         ker_lm(:,:,:,:) = 0.d0
!         ! Fast-averaging
!         allocate(PS_1_fast_lm(nfeat_1,degen,nmol_1),PS_2_fast_lm(nfeat_2,degen,nmol_2))
!         PS_1_fast_lm(:,:,:) = 0.d0
!         PS_2_fast_lm(:,:,:) = 0.d0
!         do i=1,nmol_1
!          do ii=1,int(natoms_1(i))
!           PS_1_fast_lm(:,:,i) = PS_1_fast_lm(:,:,i) + PS_1_lm(:,:,ii,i) / natoms_1(i)
!          enddo
!         enddo
!         do i=1,nmol_2
!          do ii=1,int(natoms_2(i))
!           PS_2_fast_lm(:,:,i) = PS_2_fast_lm(:,:,i) + PS_2_lm(:,:,ii,i) / natoms_2(i)
!          enddo
!         enddo
!         do i=1,nmol_1
!          j0=1
!          if (hermiticity) j0=i
!          !$OMP PARALLEL DO SHARED(ker_lm,PS_1_fast_lm,PS_2_fast_lm,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu)
!          do j=j0,nmol_2
!           do mu=1,degen
!            do nu=1,degen
!             ker_lm(i,j,mu,nu) = (dot_product(PS_1_fast_lm(:,mu,i),PS_2_fast_lm(:,nu,j)))
!             if (hermiticity) ker_lm(j,i,nu,mu) = ker_lm(i,j,mu,nu)
!            enddo
!           enddo
!          enddo
!          !$OMP END PARALLEL DO
!         enddo
!        else
!         ! Read in scalar power spectrum file(s)
!         if (psname0(1).eq.'' .or. psname0(2).eq.'') stop 'ERROR: scalar power spectrum file(s) required!'
!         do i=1,2
!          open(unit=41,file=psname0(i),status='old',access='stream',form='unformatted')
!          inquire(unit=41,size=bytes)
!          reals = bytes/8
!          if (i.eq.1) then
!           allocate(raw_PS0_1(reals))
!           read(41,pos=1) raw_PS0_1
!           close(41)
!           nfeat0_1 = reals / (nmol_1*natmax_1)
!           call C_F_POINTER(C_LOC(raw_PS0_1),PS0_1,[nfeat0_1,natmax_1,nmol_1])
!          else
!           allocate(raw_PS0_2(reals))
!           read(41,pos=1) raw_PS0_2
!           close(41)
!           nfeat0_2 = reals / (nmol_2*natmax_2)
!           call C_F_POINTER(C_LOC(raw_PS0_2),PS0_2,[nfeat0_2,natmax_2,nmol_2])
!          endif
!         enddo
!         ! Build kernel
!         allocate(ker_lm(nmol_1,nmol_2,degen,degen),ker(nmol_1,nmol_2))
!         ker_lm(:,:,:,:) = 0.d0
!         ker(:,:) = 0.d0
!         do i=1,nmol_1
!          j0=1
!          if (hermiticity) j0=i
!          !$OMP PARALLEL DO SHARED(ker_lm,ker,PS_1_lm,PS_2_lm,PS0_1,PS0_2,natoms_1,natoms_2) PRIVATE(ii,jj,mu,nu,k0)
!          do j=j0,nmol_2
!           do ii=1,int(natoms_1(i))
!            do jj=1,int(natoms_2(j))
!             k0 = (dot_product(PS0_1(:,ii,i),PS0_2(:,jj,j)))
!             do mu=1,degen
!              do nu=1,degen
!               ker_lm(i,j,mu,nu) = ker_lm(i,j,mu,nu) + dot_product(PS_1_lm(:,mu,ii,i),PS_2_lm(:,nu,jj,j)) * k0**(zeta-1)
!              enddo
!             enddo
!            enddo
!           enddo
!           ker_lm(i,j,:,:) = ker_lm(i,j,:,:) / (natoms_1(i)*natoms_2(j))
!           if (hermiticity) then
!            do mu=1,degen
!             do nu=1,degen
!              ker_lm(j,i,nu,mu) = ker_lm(i,j,mu,nu)
!             enddo
!            enddo
!           endif
!          enddo
!          !$OMP END PARALLEL DO
!         enddo
!        endif

!        deallocate(arg,natoms_1,natoms_2,raw_PS_1,raw_PS_2)
!        if (allocated(PS_1_fast_lm)) deallocate(PS_1_fast_lm)
!        if (allocated(PS_2_fast_lm)) deallocate(PS_2_fast_lm)

    endif

    ! Print kernel
    open(42,file=ofile,access='stream',form='unformatted')
     if (lm.eq.0) then
      write(42,pos=1) ker
     else
      write(42,pos=1) ker_lm
     endif
    close(42)

    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)

end program
