program sagpr_apply
    use sagpr
    use apply
    USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
    implicit none

    character(len=100), allocatable :: arg(:),keys(:)
    integer i,j,k,l,ios,numkeys,port
    character(len=100) ofile,model,fname,sock_arg(3),hostname
    integer t1,t2,cr,ts,tf
    real*8 rate
    logical readnext,use_socket

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
    numkeys = 5
    allocate(keys(numkeys))
    model = ''
    ofile = 'prediction.out'
    fname = ''
    use_socket = .false.
    keys = (/'-m  ','-o  ','-f  ','-s  ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') read(arg(i+1),'(A)') model
     if (arg(i).eq.'-o') read(arg(i+1),'(A)') ofile
     if (arg(i).eq.'-f') read(arg(i+1),'(A)') fname
     if (arg(i).eq.'-s') then
      use_socket = .true.
      readnext = .true.
      l = 0
      do k=i+1,nargs
       do j=1,numkeys
        readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keys(j))))
       enddo
       if (readnext) then
        l = l + 1
        read(arg(k),*) sock_arg(l)
       endif
      enddo
     endif
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    if (model.eq.'') stop 'ERROR: model file required!'
    if (fname.eq.'' .and. (.not. use_socket)) stop 'ERROR: file name required!'

!************************************************************************************
! GET MODEL
!************************************************************************************

    call get_model(model)

!************************************************************************************
! READ IN DATA AND MAKE PREDICTIONS
!************************************************************************************

    if (.not. use_socket) then
     ! Open xyz file
     open(unit=15,file=trim(adjustl(fname)),status="old",access="stream",form="formatted")
    else
     ! Set up sockets
     if (l.eq.0) then
     else if (l.eq.2) then
      hostname = trim(adjustl(sock_arg(1)))//achar(0)
      read(sock_arg(2),*) port
     else
      stop 'ERROR: wrong number of arguments given for socket!'
     endif
    endif

    ! Initialize the system clock
    call system_clock(count_rate=cr)
    rate = real(cr)

    ios = 1
    ! Open output file
    if (.not. use_socket) open(unit=33,file=ofile,access='stream',form='formatted')
    ! Keep going through input
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      if (.not. use_socket) then
       ! Read the next frame from xyz file
       call read_frame(15,periodic)
      else
       ! Read the next frame from socket
      endif

      call system_clock(t1)

      ! Do prediction
      call predict_frame(rate)

      ! Print predictions
      if (.not. use_socket) then
       write(33,*) (prediction_lm(1,j),j=1,degen)
       flush(33)
      else
      endif

     endif

     call system_clock(t2)
     write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
     write(*,*)
    enddo
    close(33)

    ! Array deallocation
    if (allocated(xyz)) deallocate(xyz,atname,natoms,comment,sparsification,cell,PS_tr_lam,wt,arg,keys,prediction_lm,PS)
    if (allocated(natoms_tr)) deallocate(natoms_tr)
    if (allocated(PS0)) deallocate(PS0)
    if (allocated(ker)) deallocate(ker)
    if (allocated(ker_lm)) deallocate(ker_lm)
    if (allocated(sparsification0)) deallocate(sparsification0,PS_tr_0)
    if (allocated(components)) deallocate(components)
    if (allocated(w3j)) deallocate(w3j)

end program
