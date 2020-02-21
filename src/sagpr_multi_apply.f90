program sagpr_apply
    use sagpr
    use apply
    implicit none

    character(len=100), allocatable :: arg(:),keys1(:)
    integer i,j,ios,un,u2
    character(len=100) ofile,model,fifo_in,fifo_out
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
! GET MODEL
!************************************************************************************

    call get_model(model)

!************************************************************************************
! READ IN DATA SEVERAL TIMES AND MAKE PREDICTIONS
!************************************************************************************

    fifo_in = trim(adjustl(model))//'.in'
    fifo_out = trim(adjustl(model))//'.out'
    call execute_command_line('if [ -p '//trim(adjustl(fifo_in))//' ];then rm '//trim(adjustl(fifo_in))//';fi')
    call execute_command_line('mkfifo '//trim(adjustl(fifo_in))//' -m777')
    call execute_command_line('if [ -p '//trim(adjustl(fifo_out))//' ];then rm '//trim(adjustl(fifo_out))//';fi')
    call execute_command_line('mkfifo '//trim(adjustl(fifo_out))//' -m777')

    open(newunit=un,file=trim(adjustl(fifo_in)),status="old",access="stream",form="formatted")
    open(newunit=u2,file=trim(adjustl(fifo_out)),status="old",access="stream",form="formatted")

    ! Initialize the system clock
    call system_clock(count_rate=cr)
    rate = real(cr)

    write(*,*) 'Ready to accept input at '//trim(adjustl(fifo_in))

    ios = 1
    open(unit=33,file=ofile)
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      call read_fifo(un,periodic,model)

      call system_clock(t1)

      ! Do prediction
      call predict_frame(rate)

      ! Print predictions

      do i=1,nframes
       write(33,*) (prediction_lm(i,j),j=1,degen)
       flush(33)
       write(u2,*) (prediction_lm(i,j),j=1,degen)
      enddo

     endif

     call system_clock(t2)
     write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
     write(*,*)
    enddo
    close(33)
    call execute_command_line('rm '//trim(adjustl(model))//'.{in,out}')

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
