program sagpr_apply
    use sagpr
    use apply
    implicit none

    character(len=100), allocatable :: arg(:),keys1(:)
    integer i,j,ios,un,u2
    character(len=100) ofile,model
    real*8 t1,t2

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
    allocate(keys1(3))
    model = ''
    ofile = 'prediction.out'
    keys1 = (/'-m  ','-o  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') read(arg(i+1),'(A)') model
     if (arg(i).eq.'-o') read(arg(i+1),'(A)') ofile
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    if (model.eq.'') stop 'ERROR: model file required!'

!************************************************************************************
! GET INFORMATION ABOUT THE MODEL
!************************************************************************************

    call set_defaults()
    arg = process_hyperparameters(model)

!************************************************************************************
! GET WEIGHTS FROM FILE
!************************************************************************************

    call get_model(model)

!************************************************************************************
! READ IN DATA SEVERAL TIMES AND MAKE PREDICTIONS
!************************************************************************************

    call execute_command_line('if [ -p my_fifo_in ];then rm my_fifo_in;fi')
    call execute_command_line('mkfifo my_fifo_in')
    call execute_command_line('if [ -p my_fifo_out ];then rm my_fifo_out;fi')
    call execute_command_line('mkfifo my_fifo_out')

    open(newunit=un,file="my_fifo_in",access="stream",form="formatted")
    open(newunit=u2,file="my_fifo_out",access="stream",form="formatted")

    ios = 1
    open(unit=33,file=ofile)
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      call read_fifo(un,periodic)

      call cpu_time(t1)

      ! Get power spectrum
      if (.not.do_scalar) then
       call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut,sparsification,rs,periodic,.true.)
      else
       call do_power_spectrum(xyz,atname,natoms,cell,nframes,natmax,lm,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut,sparsification,rs,periodic,.true.)
       call do_power_spectrum_scalar(xyz,atname,natoms,cell,nframes,natmax,0,nmax,lmax,rcut,sg,all_centres,all_species, &
     &     ncut0,sparsification0,rs,periodic,.true.)
      endif

      ! Get kernel
      if (allocated(natoms_tr)) deallocate(natoms_tr)
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
      if (allocated(prediction_lm)) deallocate(prediction_lm)
      allocate(prediction_lm(nmol,degen))
      if (lm.eq.0) then
       prediction_lm = do_prediction(ker,wt,meanval,degen,nframes,nmol)
      else
       prediction_lm = do_prediction(ker_lm,wt,meanval,degen,nframes,nmol)
      endif

      ! Print predictions

      do i=1,nframes
       write(33,*) (prediction_lm(i,j),j=1,degen)
       flush(33)
       write(u2,*) (prediction_lm(i,j),j=1,degen)
      enddo

     endif

     call cpu_time(t2)
     write(*,*) 'Time taken (seconds):',t2-t1
    enddo
    close(33)
    call execute_command_line('rm my_fifo')

    ! Array deallocation
    if (allocated(xyz)) deallocate(xyz,atname,natoms,comment,sparsification,cell,PS_tr_lam,wt,arg,keys1,prediction_lm,PS)
    if (allocated(natoms_tr)) deallocate(natoms_tr)
    if (allocated(PS0)) deallocate(PS0)
    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)
    if (allocated(sparsification0)) deallocate(sparsification0,PS_tr_0)
    if (allocated(components)) deallocate(components)
    if (allocated(w3j)) deallocate(w3j)

end program
