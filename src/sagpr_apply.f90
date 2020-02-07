program sagpr_apply
    use sagpr
    use apply
    implicit none

    character(len=100), allocatable :: arg(:),keys1(:)
    integer i,j
    character(len=100) ofile,fname,model

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

!************************************************************************************
! GET MODEL
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
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
    else
     call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg, &
     &     ncut,sparsification,rs,periodic,.true.)
     call do_power_spectrum_scalar(xyz,atname,natoms,cell,nframes,natmax,0,nmax0,lmax0,rcut0,sg0, &
     &     ncut0,sparsification0,rs0,periodic,.true.)
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
     if (lm.eq.0) stop 'ERROR: you have asked for a spherical power spectrum but provided lambda=0!'
     ker_lm = do_nonlinear_spherical_kernel(real(PS),PS_tr_lam,real(PS0),PS_tr_0,nframes,nmol,nfeat,nfeat, &
     &     nfeat0,nfeat0,natmax,1,dfloat(natoms),natoms_tr,zeta,.false.,degen)
    endif

    ! Get predictions
    allocate(prediction_lm(nmol,degen))
    if (lm.eq.0) then
     prediction_lm = do_prediction(ker,wt,meanval,degen,nframes,nmol)
    else
     prediction_lm = do_prediction(ker_lm,wt,meanval,degen,nframes,nmol)
    endif

    ! Print predictions
    open(unit=33,file=ofile)
    do i=1,nframes
     write(33,*) (prediction_lm(i,j),j=1,degen)
    enddo
    close(33)
    ! Array deallocation
    if (allocated(xyz)) deallocate(xyz,atname,natoms,comment,sparsification,cell,PS_tr_lam,wt,keys1,prediction_lm,PS)
    if (allocated(natoms_tr)) deallocate(natoms_tr)
    if (allocated(PS0)) deallocate(PS0)
    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)
    if (allocated(sparsification0)) deallocate(sparsification0,PS_tr_0)
    if (allocated(components)) deallocate(components)
    if (allocated(w3j)) deallocate(w3j)

end program
