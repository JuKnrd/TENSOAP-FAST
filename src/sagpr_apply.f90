program sagpr_apply
    use sagpr
    use apply
    implicit none

    character(len=100), allocatable :: arg(:),keys1(:)
    integer i,j,ios,un,u2
    character(len=100) ofile,model,fname
    integer t1,t2,cr,ts,tf
    real*8 rate

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
    model = ''
    ofile = 'prediction.out'
    fname = ''
    keys1 = (/'-m  ','-o  ','-f  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') read(arg(i+1),'(A)') model
     if (arg(i).eq.'-o') read(arg(i+1),'(A)') ofile
     if (arg(i).eq.'-f') read(arg(i+1),'(A)') fname
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    if (model.eq.'') stop 'ERROR: model file required!'
    if (fname.eq.'') stop 'ERROR: file name required!'

!************************************************************************************
! GET MODEL
!************************************************************************************

    call get_model(model)

!************************************************************************************
! READ IN DATA AND MAKE PREDICTIONS
!************************************************************************************

    open(unit=un,file=trim(adjustl(fname)),status="old",access="stream",form="formatted")

    ! Initialize the system clock
    call system_clock(count_rate=cr)
    rate = real(cr)

    ios = 1
    open(unit=33,file=ofile)
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      call read_frame(un,periodic)

      call system_clock(t1)

      ! Do prediction
      call predict_frame(rate)

      ! Print predictions
      do i=1,nframes
       write(33,*) (prediction_lm(i,j),j=1,degen)
       flush(33)
      enddo

     endif

     call system_clock(t2)
     write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
     write(*,*)
    enddo
    close(33)

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
