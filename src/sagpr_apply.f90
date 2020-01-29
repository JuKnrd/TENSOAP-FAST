program sagpr_apply
    use sagpr
    use apply
    implicit none

    integer numkeys
    character(len=100), allocatable :: arg(:),keylist(:),keys1(:),keys2(:)
    integer i,j,k,ii
    character(len=100) ofile,fname,model
    logical readnext
    character(len=3) symbol

!************************************************************************************
! GET COMMAND-LINE ARGUMENTS
!************************************************************************************

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    allocate(keys1(4))
    fname = ''
    model = ''
    ofile = 'prediction.out'
    keys1 = (/'-f  ','-m  ','-o  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-f') read(arg(i+1),'(A)') fname
     if (arg(i).eq.'-m') read(arg(i+1),'(A)') model
     if (arg(i).eq.'-o') read(arg(i+1),'(A)') ofile
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    if (fname.eq.'') stop 'ERROR: filename required!'
    if (model.eq.'') stop 'ERROR: model file required!'

    ! Read in hyperparameters
    call set_defaults()
    arg = process_hyperparameters(model)

!************************************************************************************
! GET INFORMATION ABOUT THE MODEL
!************************************************************************************

!    ! Get hyperparameters from model file
!    allocate(keys2(11))
!    keys2 = (/'-lm ','-n  ','-l  ','-rc ','-sg ','-c  ','-s  ','-rs ','-p  ','-z  ','NULL'/)
!    do i=1,nargs
!     arg(i) = trim(adjustl(arg(i)))
!     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
!     if (arg(i).eq.'-n') read(arg(i+1),*) nmax
!     if (arg(i).eq.'-l') read(arg(i+1),*) lmax
!     if (arg(i).eq.'-rc') read(arg(i+1),*) rcut
!     if (arg(i).eq.'-sg') read(arg(i+1),*) sg
!     if (arg(i).eq.'-z') read(arg(i+1),*) zeta
!     if (arg(i).eq.'-c') then
!      readnext = .true.
!      do k=i+1,nargs
!       do j=1,numkeys
!        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keylist(j))))
!       enddo
!       if (readnext) then
!        read(arg(k),*) symbol
!        do ii=1,nelements
!         if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_centres(ii) = .true.
!        enddo
!       endif
!      enddo
!     endif
!     if (arg(i).eq.'-s') then
!      readnext = .true.
!      do k=i+1,nargs
!       do j=1,numkeys
!        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keylist(j))))
!       enddo
!       if (readnext) then
!        read(arg(k),*) symbol
!        do ii=1,nelements
!         if (trim(adjustl(atomic_names(ii))).eq.trim(adjustl(symbol))) all_species(ii) = .true.
!        enddo
!       endif
!      enddo
!     endif
!     if (arg(i).eq.'-rs') then
!      read(arg(i+1),*) rs(1)
!      read(arg(i+2),*) rs(2)
!      read(arg(i+3),*) rs(3)
!     endif
!     if (arg(i).eq.'-p') periodic=.true.
!    enddo
!    ! Set degeneracy
!    degen = 2*lm + 1

!************************************************************************************
! GET WEIGHTS FROM FILE
!************************************************************************************

    call get_model(model)

!************************************************************************************
! READ IN DATA
!************************************************************************************

    call read_xyz(fname,periodic)

!************************************************************************************
! USE THE MODEL AND THE DATA TO MAKE A PREDICTION
!************************************************************************************

    ! Get power spectrum
    if (.not.do_scalar) then
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut,sparsification,rs,periodic,.true.)
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
    deallocate(xyz,atname,natoms,comment,sparsification,cell,PS,PS_tr_lam,natoms_tr,wt,arg)
    if (allocated(PS_tr_0)) deallocate(PS_tr_0)
    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)

end program
