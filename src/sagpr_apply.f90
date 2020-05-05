program sagpr_apply
    use sagpr
    use apply
    USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
    implicit none

    character(len=100), allocatable :: arg(:),keys(:)
    integer i,j,k,l,ios,numkeys,port,socket,inet,cbuf
    integer nat
    character(len=100) ofile,model,fname,sock_arg(3)
    character(len=1024) hostname
    character(len=12) header
    character(len=2048) initbuffer
    integer, parameter :: MSGLEN=12
    integer t1,t2,cr,ts,tf
    real*8 rate,mtxbuf(9),cell_h(3,3)
    real*8, allocatable :: msgbuffer(:)
    real*8, parameter :: bohrtoangstrom = 0.52917721d0
    logical readnext,use_socket,hasdata

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
    numkeys = 6
    allocate(keys(numkeys))
    model = ''
    ofile = 'prediction.out'
    fname = ''
    use_socket = .false.
    inet = 1
    keys = (/'-m  ','-o  ','-f  ','-s  ','-u  ','NULL'/)
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
     if (arg(i).eq.'-u') inet = 0
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

    ! Open xyz file
    open(unit=15,file=trim(adjustl(fname)),status="old",access="stream",form="formatted")
    if (use_socket) then
     ! Get atom names
     call read_frame(15,periodic)
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
    ! Open output file
    if (.not. use_socket) open(unit=33,file=ofile,access='stream',form='formatted')
    ! Keep going through input
    do while (ios.ne.0)
     open(unit=73,file='EXIT',status='old',iostat=ios)
     if (ios.ne. 0) then

      if (.not. use_socket) then
       ! Read the next frame from xyz file
       call read_frame(15,periodic)

       call system_clock(t1)

       ! Do prediction
       call predict_frame(rate)

       ! Print predictions
       write(33,*) (prediction_lm(1,j),j=1,degen)
       flush(33)

       call system_clock(t2)
       write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
       write(*,*)

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
        if (allocated(xyz)) deallocate(xyz,natoms,cell)
        call readbuffer(socket, mtxbuf, 9)  ! Cell matrix
        cell_h = reshape(mtxbuf, (/3,3/))
        cell_h = transpose(cell_h)
        call readbuffer(socket, mtxbuf, 9)
        call readbuffer(socket, cbuf)
        nat = cbuf
        natmax = nat
        nframes = 1
        allocate(xyz(nframes,natmax,3),natoms(nframes),msgbuffer(3*nat),cell(nframes,3,3))
        natoms(1) = nat
        call readbuffer(socket, msgbuffer, nat*3) 
        ! Read in atoms
        do i=1,nat
         xyz(1,i,:) = msgbuffer(3*(i-1)+1:3*i) * bohrtoangstrom
        enddo
        cell(1,:,:) = cell_h(:,:) * bohrtoangstrom

        ! With all of the necessary data read in, do the prediction
        call system_clock(t1)

        ! Do prediction
        call predict_frame(rate)
        write(*,*) prediction_lm

        call system_clock(t2)
        write(*,'(A,F6.3,A)') '===>Time taken: ',(t2-t1)/rate,' seconds'
        write(*,*)

        ! We have the data
        hasdata = .true.
        deallocate(xyz,atname,natoms,msgbuffer,cell)
       elseif (trim(header)=='GETFORCE') then
        write(*,*) 'EXPECTING FORCE'
        stop
       else
        write(*,*) 'ERROR: unexpected header ',header
        stop
       endif
      endif
     endif
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
