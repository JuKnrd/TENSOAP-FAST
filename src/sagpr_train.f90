program sagpr_train
    use iso_c_binding
    implicit none

    integer nargs,numkeys,ios
    character(len=100), allocatable :: arg(:),keylist(:)

    integer lm,sel(2),rdm,i,j,nmol,reals,bytes,degen,nenv,ntr,mu,nu,nte,k
    real*8 rg,ftr,mean_val
    character(len=100) feat,ker(2),weights
    logical peratom,predict,sparsify,readnext
    real*8, allocatable, target :: raw_K1(:), raw_K2(:)
    real*8, pointer :: K1(:,:),K2(:,:),K1_lm(:,:,:,:),K2_lm(:,:,:,:)
    real*8, allocatable :: data_points(:,:),ktr(:,:),data_vec(:)
    real*8, allocatable :: ktr_NM(:,:),ktr_MM(:,:),ktr_MNNM(:,:),data_vec_MN(:),wt(:),kte(:,:),pred(:),final_pred(:,:)
    integer, allocatable :: train_list(:),all_list(:),test_list(:)

    integer, parameter :: lwmax = 10000
    integer info,lwork
    real*8 work(lwmax)
    real*8, allocatable :: S(:),U(:,:),VT(:,:),Smat(:,:)
    real*8, parameter :: threshold=1.d-15

    ! Get input arguments
    nargs = iargc()
    allocate(arg(nargs+1))
    do i=1,nargs
     call getarg(i,arg(i))
    enddo
    arg(nargs+1) = 'NULL'

    ! Parse these arguments
    lm = 0
    rg = -1.d-10
    ftr = 1.d0
    feat = ''
    ker = (/'',''/)
    sparsify = .false.
    sel = (/0,0/)
    rdm = 0
    weights = ''
    peratom = .false.
    predict = .false.
    numkeys = 12
    allocate(keylist(numkeys))
    keylist=(/'-lm ','-rg ','-ft ','-f  ', '-k  ','-sf ','-sel','-rdm','-w  ','-pat','-pr ','NULL'/)
    do i=1,nargs
     arg(i) = trim(adjustl(arg(i)))
     if (arg(i).eq.'-lm') read(arg(i+1),*) lm
     if (arg(i).eq.'-rg') read(arg(i+1),*) rg
     if (arg(i).eq.'-ft') read(arg(i+1),*) ftr
     if (arg(i).eq.'-f')  read(arg(i+1),*) feat
     if (arg(i).eq.'-k') then
      read(arg(i+1),*) ker(1)
      readnext = .true.
      do j=1,numkeys
       readnext = readnext.and.(arg(i+2).ne.trim(adjustl(keylist(j))))
      enddo
      if (readnext) read(arg(i+2),*) ker(2)
     endif
     if (arg(i).eq.'-sf') sparsify=.true.
     if (arg(i).eq.'-sel') then
      read(arg(i+1),*) sel(1)
      readnext = .true.
      do j=1,numkeys
       readnext = readnext.and.(arg(i+2).ne.trim(adjustl(keylist(j))))
      enddo
      if (readnext) read(arg(i+2),*) sel(2)
     endif
     if (arg(i).eq.'-rdm') read(arg(i+1),*) rdm
     if (arg(i).eq.'-w') read(arg(i+1),*) weights
     if (arg(i).eq.'-pat') peratom=.true.
     if (arg(i).eq.'-pr') predict=.true.
    enddo

    ! Check for arguments that are required
    if (rg.lt.0.d0) stop 'ERROR: regularization required!'
    if (feat.eq.'') stop 'ERROR: features file required!'
    if (ker(1).eq.'') stop 'ERROR: kernel required!'
    if (sparsify.and.(ker(2).eq.'')) stop 'ERROR: two kernels required!'
    if ((sel(1).eq.0).and.(sel(2).eq.0).and.(rdm.eq.0)) stop 'ERROR: either a selection or randomization required!'
    if ((rdm.eq.0).and.(sel(2).le.sel(1))) stop 'ERROR: selection given incorrectly!'
    if (weights.eq.'') stop 'ERROR: weights file name required!'

    ! Read in features
    nmol = 0
    degen = 2*lm + 1
    open(unit=30,file=feat,status='old')
    do
     read(30,*,iostat=ios)
     if (ios.ne.0) exit
     nmol = nmol + 1
    enddo
    close(30)
    write(*,*) 'Found',nmol,'data points'
    allocate(data_points(nmol,degen))
    open(unit=30,file=feat,status='old')
    do i=1,nmol
     read(30,*) (data_points(i,j),j=1,degen)
    enddo

    ! Read in kernel or kernels
    open(unit=31,file=ker(1),status='old',access='stream',form='unformatted')
    inquire(unit=31,size=bytes)
    reals = bytes/8
    allocate(raw_K1(reals))
    read(31,pos=1) raw_K1
    close(31)
    if (lm.eq.0) then
     nenv = reals / nmol
     call C_F_POINTER(C_LOC(raw_K1),K1,[nmol,nenv])
    else
     nenv = reals / (nmol * degen * degen)
     call C_F_POINTER(C_LOC(raw_K1),K1_lm,[nmol,nenv,degen,degen])
    endif
    write(*,*) 'Read in kernel',ker(1)
    if (ker(2).ne.''.and.sparsify) then
     open(unit=31,file=ker(2),status='old',access='stream',form='unformatted')
     inquire(unit=31,size=bytes)
     reals = bytes/8
     allocate(raw_K2(reals))
     read(31,pos=1) raw_K2
     close(31)
     write(*,*) 'Read in kernel',ker(2)
     if (lm.eq.0) then
      call C_F_POINTER(C_LOC(raw_K2),K2,[nenv,nenv])
     else
      call C_F_POINTER(C_LOC(raw_K2),K2_lm,[nenv,nenv,degen,degen])
     endif
    endif

    ! Set up the matrix to be inverted
    if (rdm.ne.0) then
     ntr = rdm
     allocate(train_list(ntr),all_list(nmol))
     do i=1,nmol
      all_list(i) = i
     enddo
     call shuffle(all_list,nmol)
     train_list = all_list(1:ntr)
    else
     ntr = 1 + sel(2) - sel(1)
     allocate(train_list(ntr))
     do i=1,ntr
      train_list(i) = sel(1) + i - 1
     enddo
    endif

    if (.not.sparsify) then
     if (lm.eq.0) then
      allocate(ktr(ntr,ntr))
      do i=1,ntr
       do j=1,ntr
        ktr(i,j) = K1(train_list(i),train_list(j))
       enddo
       ktr(i,i) = ktr(i,i) + rg
      enddo
     else
      allocate(ktr(ntr*degen,ntr*degen))
      do i=1,ntr
       do mu=1,degen
        do j=1,ntr
         do nu=1,degen
          ktr((i-1)*degen + mu,(j-1)*degen + nu) = K1_lm(train_list(i),train_list(j),mu,nu)
         enddo
        enddo
        ktr((i-1)*degen + mu,(i-1)*degen + mu) = ktr((i-1)*degen + mu,(i-1)*degen + mu) + rg
       enddo
      enddo
     endif
    else
     ! Deal with sparsification kernels
     if (lm.eq.0) then
      allocate(ktr_NM(ntr,nenv),ktr_MM(nenv,nenv))
      do i=1,ntr
       do j=1,nenv
        ktr_NM(i,j) = K1(train_list(i),j)
       enddo
      enddo
      do i=1,nenv
       do j=1,nenv
        ktr_MM(i,j) = K2(i,j)
       enddo
      enddo
     else
       allocate(ktr_NM(ntr*degen,nenv*degen),ktr_MM(nenv*degen,nenv*degen))
       do i=1,ntr
        do mu=1,degen
         do j=1,nenv
          do nu=1,degen
           ktr_NM((i-1)*degen + mu,(j-1)*degen + nu) = K1_lm(train_list(i),j,mu,nu)
          enddo
         enddo
        enddo
       enddo
       do i=1,nenv
        do mu=1,degen
         do j=1,nenv
          do nu=1,degen
           ktr_MM((i-1)*degen + mu,(j-1)*degen + nu) = K2_lm(i,j,mu,nu)
          enddo
         enddo
        enddo
       enddo
     endif
     allocate(ktr_MNNM(nenv*degen,nenv*degen))
     ktr_MNNM = matmul(transpose(ktr_NM),ktr_NM) + rg*ktr_MM
    endif

    ! Reshape and/or transform data
    mean_val = 0.d0
    if (lm.eq.0) then
     allocate(data_vec(ntr))
     do i=1,ntr
      data_vec(i) = data_points(i,1)
     enddo
     do i=1,ntr
      mean_val = mean_val + (data_vec(i) / float(ntr))
     enddo
     do i=1,ntr
      data_vec(i) = data_vec(i) - mean_val
     enddo
    else
     allocate(data_vec(ntr*degen))
     do i=1,ntr
      do mu=1,degen
       data_vec((i-1)*degen + mu) = data_points(i,mu)
      enddo
     enddo
    endif
    if (sparsify) then
     allocate(data_vec_MN(nenv*degen))
     data_vec_MN = matmul(transpose(ktr_NM),data_vec)
    endif

    ! Solve linear problem
    write(*,*) 'Doing regression'
    if (.not.sparsify) then
     allocate(S(ntr),U(ntr,ntr),VT(ntr,ntr))
     lwork = -1
     call DGESVD('All','All',ntr,ntr,ktr,ntr,S,U,ntr,VT,ntr,work,lwork,info)
     lwork = min(lwmax,int(work(1)))
     call DGESVD('All','All',ntr,ntr,ktr,ntr,S,U,ntr,VT,ntr,work,lwork,info)
     allocate(Smat(ntr,ntr))
     Smat(:,:) = 0.d0
     do i=1,ntr
      if (S(i).ge.threshold) Smat(i,i) = 1.d0 / S(i)
     enddo
     ktr = matmul(transpose(VT),matmul(Smat,transpose(U)))
    else
     allocate(S(nenv),U(nenv,nenv),VT(nenv,nenv))
     lwork = -1
     call DGESVD('All','All',nenv,nenv,ktr_MNNM,nenv,S,U,ntr,VT,nenv,work,lwork,info)
     lwork = min(lwmax,int(work(1)))
     call DGESVD('All','All',nenv,nenv,ktr_MNNM,nenv,S,U,ntr,VT,nenv,work,lwork,info)
     allocate(Smat(ntr,ntr))
     Smat(:,:) = 0.d0
     do i=1,ntr
      if (S(i).ge.threshold) Smat(i,i) = 1.d0 / S(i)
     enddo
     ktr_MNNM = matmul(transpose(VT),matmul(Smat,transpose(U)))
    endif

    ! Print weights
    if (.not.sparsify) then
     allocate(wt(0:ntr))
     wt(0) = mean_val
     wt(1:ntr) = matmul(ktr,data_vec)
    else
     allocate(wt(0:nenv))
     wt(0) = mean_val
     wt(1:nenv) = matmul(ktr_MNNM,data_vec)
    endif
    open(32,file=weights,access='stream',form='unformatted')
    write(32,pos=1) wt
    close(32)

    ! Print training set
    open(33,file='training_set_'//trim(adjustl(weights))//'.txt')
    do i=1,ntr
     write(33,*) train_list(i)
    enddo
    close(33)

    if (predict) then
     write(*,*) 'Doing predictions'
     nte = nmol-ntr
     allocate(test_list(nte))
     j=0
     do i=1,nmol
      if (.not.any(train_list==i)) then
       j = j + 1
       test_list(j) = i
      endif
     enddo
     if (.not.sparsify) then
      if (lm.eq.0) then
       allocate(kte(nte,ntr))
       do i=1,nte
        do j=1,ntr
         kte(i,j) = K1(test_list(i),train_list(j))
        enddo
       enddo
      else
       allocate(kte(nte*degen,ntr*degen))
       do i=1,nte
        do mu=1,degen
         do j=1,ntr
          do nu=1,degen
           kte((i-1)*degen + mu,(j-1)*degen + nu) = K1_lm(test_list(i),train_list(j),mu,nu)
          enddo
         enddo
        enddo
       enddo
      endif
     else
      if (lm.eq.0) then
       allocate(kte(nte,nenv))
       do i=1,nte
        do j=1,nenv
         kte(i,j) = K1(train_list(i),j)
        enddo
       enddo
      else
       allocate(kte(nte*degen,nenv*degen))
       do i=1,nte
        do mu=1,degen
         do j=1,nenv
          do nu=1,degen
           kte((i-1)*degen + mu,(j-1)*degen + nu) = K1_lm(train_list(i),j,mu,nu)
          enddo
         enddo
        enddo
       enddo
      endif
     endif
     allocate(pred(nte*degen))
     pred = matmul(kte,wt(1:size(wt))) + mean_val
     ! Reshape the array if needed
     allocate(final_pred(nte,degen))
     k = 0
     do i=1,nte
      do j=1,degen
       k = k + 1
       final_pred(i,j) = pred(k)
      enddo
     enddo
     open(unit=34,file='predictions_'//trim(adjustl(weights))//'.txt')
     do i=1,nte
      write(34,*) (final_pred(i,j),j=1,degen)
     enddo
     close(34)
     deallocate(test_list,kte,pred,final_pred)
    endif

    !TODO: deallocate
    deallocate(S,U,VT,Smat,train_list,wt,data_points)
    if (allocated(raw_K1)) deallocate(raw_K1)
    if (allocated(raw_K2)) deallocate(raw_K2)
    if (allocated(ktr)) deallocate(ktr)
    if (allocated(data_vec)) deallocate(data_vec)
    if (allocated(ktr_NM)) deallocate(ktr_NM)
    if (allocated(ktr_MM)) deallocate(ktr_MM)
    if (allocated(ktr_MNNM)) deallocate(ktr_MNNM)
    if (allocated(data_vec_MN)) deallocate(data_vec_MN)

end program

subroutine shuffle(list,ln)
    implicit none

    integer ln
    integer, intent(inout) :: list(ln)
    integer :: i,randpos,temp
    real :: r

    do i=size(list),2,-1
     call random_number(r)
     randpos = int(r*i) + 1
     temp = list(randpos)
     list(randpos) = list(i)
     list(i) = temp
    enddo

end subroutine
