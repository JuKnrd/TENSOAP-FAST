program sagpr_apply
    use apply
    use io
    USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
    implicit none

    type(SAGPR_Model), allocatable :: GPR(:)
    type(Frame_XYZ) :: frames
    character(len=100), allocatable :: arg(:),keys(:)
    integer i,j,k,l,ios,numkeys,port,socket,inet,cbuf,nargs,fr
    integer nat,s_frame,n_mod
    character(len=100) fname,sock_arg(3),ifile
    character(len=100), allocatable :: model(:),ofile(:),lodeparams(:)
    character(len=1024) hostname
    character(len=12) header
    character(len=2048) initbuffer
    integer, parameter :: MSGLEN=12
    integer t1,t2,cr,ts,tf
    real*8 rate
    logical readnext,use_socket,hasdata,verbose,atomic,exists,fixed_cell

namelist/input/model,fname,ofile,use_socket,sock_arg,inet,verbose,atomic,s_frame

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

    ! Do a first pass to find the number of input arguments
    numkeys = 12
    allocate(keys(numkeys))
    keys = (/'-m  ','-o  ','-f  ','-s  ','-u  ','-v  ','-a  ','-b  ','-i  ','-l  ','-fc ','NULL'/)
    n_mod = -1
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') then
      n_mod = 0
      readnext = .true.
      do k=i+1,nargs
       if (readnext) then
        do j=1,numkeys
         readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keys(j))))
        enddo
        if (readnext) n_mod = n_mod + 1
       endif
      enddo
     endif
    enddo
    if (n_mod.eq.-1) stop 'ERROR: no model given!'
    allocate(GPR(n_mod),model(n_mod))

    ! Check the number of output files given is also correct
    l = -1
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-o') then
      l = 0
      readnext = .true.
      do k=i+1,nargs
       if (readnext) then
        do j=1,numkeys
         readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keys(j))))
        enddo
        if (readnext) l = l + 1
       endif
      enddo
     endif
    enddo
    if (l.eq.-1) stop 'ERROR: no output files given!'
    if (l.ne.n_mod) stop 'ERROR: wrong number of outputs given!'
    allocate(ofile(n_mod))

    ! Check number of LODE inputs is either zero or equal to the number of models
    l = -1
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-l') then
      l = 0
      readnext = .true.
      do k=i+1,nargs
       if (readnext) then
        do j=1,numkeys
         readnext = readnext.and.(trim(adjustl(arg(k))).ne.trim(adjustl(keys(j))))
        enddo
        if (readnext) l = l + 1
       endif
      enddo
     endif
    enddo
    if (l.gt.-1 .and. l.ne.n_mod) stop 'ERROR: wrong number of LODE files given!'
    allocate(lodeparams(n_mod))
    if (l.eq.-1) then
     do i=1,n_mod
      lodeparams(i) = 'NONE'
     enddo
    endif

    ! Parse these arguments
    ifile = 'tensoap.in'
    model(:) = ''
    ofile(:) = 'prediction.out'
    fname = ''
    use_socket = .false.
    inet = 1
    GPR(:)%verbose = .false.
    GPR(:)%atomic = .false.
    s_frame = -1
    fixed_cell = .false.
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-m') then
      do k=1,n_mod
       read(arg(i+k),'(A)') model(k)
      enddo
     endif
     if (arg(i).eq.'-o') then
      do k=1,n_mod
       read(arg(i+k),'(A)') ofile(k)
      enddo
     endif
     if (arg(i).eq.'-l') then
      do k=1,n_mod
       read(arg(i+k),'(A)') lodeparams(k)
      enddo
     endif
     if (arg(i).eq.'-f') read(arg(i+1),'(A)') fname
     if (arg(i).eq.'-b') read(arg(i+1),'(I10)') s_frame
     if (arg(i).eq.'-i') read(arg(i+1),'(A)') ifile
     if (arg(i).eq.'-v') GPR(:)%verbose=.true.
     if (arg(i).eq.'-a') GPR(:)%atomic=.true.
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
     if (arg(i).eq.'-u') inet = 0
     if (arg(i).eq.'-fc') fixed_cell=.true.
    enddo
    deallocate(arg)

    ! Check for arguments that are required
    do k=1,n_mod
     if (model(k).eq.'') stop 'ERROR: model file required!'
    enddo
    if (fname.eq.'') stop 'ERROR: file name required!'

!************************************************************************************
! GET MODEL
!************************************************************************************

    do k=1,n_mod
     call get_model(GPR(k),model(k))
     call get_LODE(GPR(k),lodeparams(k),fixed_cell)
    enddo

    write(*,*) GPR(1)%lode_params%fixed_cell

!************************************************************************************
! READ IN DATA AND MAKE PREDICTIONS
!************************************************************************************

    ! Open xyz file
    open(unit=15,file=trim(adjustl(fname)),status="old",access="stream",form="formatted",iostat=ios)
    if (ios.ne.0) stop 'ERROR: input file does not exist!'
    if (use_socket) then
     ! Get atom names
     call read_frame(frames,15,GPR(1)%verbose,GPR(1)%periodic)
     ! Set up socket
     hostname = trim(adjustl(sock_arg(1)))//achar(0)
     read(sock_arg(2),*) port
     call open_socket(socket,inet,port,hostname)
    endif

    ! Initialize the system clock
    call system_clock(count_rate=cr)
    rate = real(cr)

    ios = 1
    hasdata = .false.
    fr = 0
    ! Open output file
    if (.not. use_socket) then
     do l=1,n_mod
      if (s_frame.le.0) then
       open(unit=32+l,file=ofile(l),access='stream',form='formatted')
      else
       inquire(file=ofile(l),exist=exists)
       if (exists) then
        open(unit=32+l,file=ofile(l),access='stream',form='formatted',status='old',position='append')
       else
        open(unit=32+l,file=ofile(l),access='stream',form='formatted')
       endif
      endif
     enddo
    endif
    if (use_socket .and. GPR(1)%atomic) then
     do l=1,n_mod
      open(unit=32+l,file=ofile(l),access='stream',form='formatted')
     enddo
    endif
    ! Keep going through input
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      if (.not. use_socket) then
       ! Read the next frame from xyz file
       call read_frame(frames,15,GPR(1)%verbose,GPR(1)%periodic)

       fr = fr + 1
       if (fr.ge.s_frame) then
         call system_clock(t1)

         do l=1,n_mod
          ! Do prediction
          call predict_frame(GPR(l),frames,rate)

          ! Rescale predictions about the mean
          call rescale_predictions(GPR(l))

          ! Print predictions
          call print_predictions(GPR(l),32+l)

          call system_clock(t2)
          if (GPR(1)%verbose) write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
          if (GPR(1)%verbose) write(*,'(A,I0)') '===> Frame ',fr
          if (GPR(1)%verbose) write(*,*)
         enddo
       endif

      else
       ! Read from the socket
       call readbuffer(socket,header,MSGLEN)
       if (trim(header)=='STATUS') then
        ! The wrapper is inquiring what we are doing
        if (hasdata) then
         call writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
        else
         call writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
        endif
       elseif (trim(header)=='INIT') then
        ! The wrapper is sending an initialization string (which we will ignore)
        call readbuffer(socket, cbuf)
        call readbuffer(socket, initbuffer, cbuf)
       elseif (trim(header)=='POSDATA') then
        ! The wrapper is sending atom positions
        call read_frame_socket(frames,nat,socket)

        ! With all of the necessary data read in, do the prediction
        call system_clock(t1)

        do l=1,n_mod
         ! Do prediction
         call predict_frame(GPR(l),frames,rate)

         ! Rescale committee predictions about their mean
         call rescale_predictions(GPR(l))
        enddo

        call system_clock(t2)
        if (GPR(1)%verbose) write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
        if (GPR(1)%verbose) write(*,*)

        ! We have the data
        hasdata = .true.
        deallocate(frames%xyz,frames%natoms,frames%cell)
       elseif (trim(header)=='GETFORCE') then
        ! Send information back to the wrapper
        call send_information_socket(GPR,nat,socket)
        hasdata = .false.
       else
        write(*,*) 'ERROR: unexpected header ',header
        stop
       endif
      endif
     endif
    enddo
    close(33)

    ! Array deallocation
    if (allocated(frames%xyz)) deallocate(frames%xyz,frames%atname, &
     &     frames%natoms,frames%comment,frames%cell)
    if (allocated(frames%xyz)) deallocate(GPR(1)%PS_tr_lam,GPR(1)%wt,arg,keys,GPR(1)%prediction_lm,GPR(1)%PS,GPR(1)%sparsification)
    if (allocated(GPR(1)%natoms_tr)) deallocate(GPR(1)%natoms_tr)
    if (allocated(GPR(1)%PS0)) deallocate(GPR(1)%PS0)
    if (allocated(GPR(1)%ker)) deallocate(GPR(1)%ker)
    if (allocated(GPR(1)%ker_lm)) deallocate(GPR(1)%ker_lm)
    if (allocated(GPR(1)%sparsification0)) deallocate(GPR(1)%sparsification0,GPR(1)%PS_tr_0)
    if (allocated(GPR(1)%components)) deallocate(GPR(1)%components)
    if (allocated(GPR(1)%w3j)) deallocate(GPR(1)%w3j)
    if (allocated(GPR(1)%meanval_c)) deallocate(GPR(1)%meanval_c,GPR(1)%wt_c,GPR(1)%prediction_lm_c)
    if (allocated(GPR(1)%atname_at)) deallocate(GPR(1)%atname_at)
    do i=1,n_mod
     deallocate(GPR,model,ofile)
    enddo

end program
