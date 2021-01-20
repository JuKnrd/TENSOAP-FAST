module io
 use apply

 contains

!****************************************************************************************************************
subroutine get_model(GPR,model)
 implicit none

  type(SAGPR_Model), intent(inout) :: GPR
  character(len=100) model
  real*8, allocatable, target :: raw_model(:)
  real*8 a1,a2
  integer i,j,k,l,ncen,nspec,reals,bytes

   ! Read in power spectrum file(s)
   open(unit=41,file=trim(adjustl(model)),status='old',access='stream',form='unformatted')
   inquire(unit=41,size=bytes)
   reals = bytes/8
   allocate(raw_model(reals))
   read(41,pos=1) raw_model
   GPR%do_scalar = (raw_model(1).ne.0.d0)
   GPR%lm = int(raw_model(2))
   GPR%degen=2*GPR%lm + 1
   GPR%zeta = int(raw_model(3))
   GPR%periodic = (int(raw_model(4)).eq.1)
   GPR%nmol = int(raw_model(5))
   GPR%nfeat = int(raw_model(6))
   i = 6
   if (GPR%do_scalar) then
    GPR%nfeat0 = int(raw_model(7))
    i = 7
   endif
   i = i + 1
   GPR%nw = int(raw_model(i))
   allocate(GPR%PS_tr_lam(GPR%nmol,1,GPR%degen,GPR%nfeat))
   do j=1,GPR%nmol
    do k=1,GPR%degen
     do l=1,GPR%nfeat
      i = i + 1
      GPR%PS_tr_lam(j,1,k,l) = raw_model(i)
     enddo
    enddo
   enddo
   if (GPR%do_scalar) then
    allocate(GPR%PS_tr_0(GPR%nmol,1,1,GPR%nfeat0))
    do j=1,GPR%nmol
     do k=1,GPR%nfeat0
      i = i + 1
      GPR%PS_tr_0(j,1,1,k) = raw_model(i)
     enddo
    enddo
   endif

   ! Get sparsification details
   i = i + 1
   GPR%ncut = int(raw_model(i))
   allocate(GPR%sparsification(2,GPR%ncut,GPR%ncut))
   do j=1,GPR%ncut
    i = i + 1
    GPR%sparsification(1,j,1) = raw_model(i)
   enddo
   do j=1,GPR%ncut
    do k=1,GPR%ncut
     i = i + 1
     a1 = raw_model(i)
     i = i + 1
     a2 = raw_model(i)
     GPR%sparsification(2,j,k) = dcmplx(a1,a2)
    enddo
   enddo
   if (GPR%do_scalar) then
    i = i + 1
    GPR%ncut0 = int(raw_model(i))
    allocate(GPR%sparsification0(2,GPR%ncut0,GPR%ncut0))
    do j=1,GPR%ncut0
     i = i + 1
     GPR%sparsification0(1,j,1) = raw_model(i)
    enddo
    do j=1,GPR%ncut0
     do k=1,GPR%ncut0
      i = i + 1
      a1 = raw_model(i)
      i = i + 1
      a2 = raw_model(i)
      GPR%sparsification0(2,j,k) = dcmplx(a1,a2)
     enddo
    enddo
   endif

   if (GPR%ncut.ne.GPR%nfeat) stop 'ERROR: ncut .ne. nfeat!'
   if (GPR%do_scalar) then
    if (GPR%ncut0.ne.GPR%nfeat0) stop 'ERROR: ncut0 .ne. nfeat0!'
   endif
   
   ! Get weights
   allocate(GPR%meanval_c(GPR%nw),GPR%wt_c(GPR%nmol*GPR%degen,GPR%nw))
   do k=1,GPR%nw
    i = i + 1
    GPR%meanval_c(k) = raw_model(i)
   enddo
   do k=1,GPR%nw
    do j=1,GPR%nmol*GPR%degen
     i = i + 1
     GPR%wt_c(j,k) = raw_model(i)
    enddo
   enddo

   ! Get hyperparameters
   i = i + 1
   GPR%nmax = int(raw_model(i))
   if (GPR%nmax.eq.-1) GPR%nmax=nmax_default
   i = i + 1
   GPR%lmax = int(raw_model(i))
   if (GPR%lmax.eq.-1) GPR%lmax=lmax_default
   i = i + 1
   GPR%rcut = raw_model(i)
   if (GPR%rcut.lt.0.d0) GPR%rcut=rcut_default
   i = i + 1
   GPR%sg = raw_model(i)
   if (GPR%sg.lt.0.d0) GPR%sg=sg_default
   i = i + 1
   GPR%rs(1) = raw_model(i)
   i = i + 1
   GPR%rs(2) = raw_model(i)
   i = i + 1
   GPR%rs(3) = raw_model(i)
   if (GPR%rs(3).eq.0.d0) GPR%rs=rs_default
   if (GPR%do_scalar) then
    i = i + 1
    GPR%nmax0 = int(raw_model(i))
    if (GPR%nmax0.eq.-1) GPR%nmax0=nmax_default
    i = i + 1
    GPR%lmax0 = int(raw_model(i))
    if (GPR%lmax0.eq.-1) GPR%lmax0=lmax_default
    i = i + 1
    GPR%rcut0 = raw_model(i)
    if (GPR%rcut0.lt.0.d0) GPR%rcut0=rcut_default
    i = i + 1
    GPR%sg0 = raw_model(i)
    if (GPR%sg0.lt.0.d0) GPR%sg0=sg_default
    i = i + 1
    GPR%rs0(1) = raw_model(i)
    i = i + 1
    GPR%rs0(2) = raw_model(i)
    i = i + 1
    GPR%rs0(3) = raw_model(i)
    if (GPR%rs0(3).eq.0.d0) GPR%rs0=rs_default
   endif
   i = i + 1
   ncen = int(raw_model(i))
   if (ncen.gt.0) then
    do j=1,ncen
     i = i + 1
     GPR%all_centres(int(raw_model(i))) = .true.
    enddo
   endif
   i = i + 1
   nspec = int(raw_model(i))
   if (nspec.gt.0) then
    do j=1,nspec
     i = i + 1
     GPR%all_species(int(raw_model(i))) = .true.
    enddo
   endif
   i = i + 1
   GPR%nu = raw_model(i)**0.5

   if (i.ne.reals) stop 'ERROR: different file size to that expected for model!'

end subroutine

!****************************************************************************************************************

subroutine read_frame(GPR,un)
 implicit none

  type(SAGPR_Model), intent(inout) :: GPR
  integer ios,un,nat
  character(len=1000) line,c1
  character(len=100) model
  integer i,j,ii

   GPR%nframes = 1
   read(un,'(A)',iostat=ios) line
   if (ios.ne.0) then
    if (GPR%verbose) write(*,*) 'End-of-file detected'
    stop
   endif
   read(line,*,iostat=ios) nat
   if (ios.ne.0) then
    call system('rm '//trim(adjustl(model))//'.{in,out}')
    write(*,*) 'Non-standard input detected; now stopping'
    stop
   endif
   GPR%natmax = nat
   if (allocated(GPR%xyz)) deallocate(GPR%xyz,GPR%atname,GPR%natoms,GPR%comment)
   allocate(GPR%xyz(GPR%nframes,GPR%natmax,3),GPR%atname(GPR%nframes,GPR%natmax), &
     &     GPR%natoms(GPR%nframes),GPR%comment(GPR%nframes))
   GPR%natoms(1) = nat
   read(un,'(A)') GPR%comment(1)
   do j=1,GPR%natoms(1)
    read(un,*) GPR%atname(1,j),(GPR%xyz(1,j,ii),ii=1,3)
   enddo

  ! Get cell data
  if (allocated(GPR%cell)) deallocate(GPR%cell)
  allocate(GPR%cell(GPR%nframes,3,3))
  if (.not.GPR%periodic) then
   GPR%cell(:,:,:) = 0.d0
  else
   do i=1,GPR%nframes
    ios = index(GPR%comment(i),'Lattice')
    if (ios.eq.0) stop 'ERROR: input file is not periodic!'
    c1 = GPR%comment(i)
    c1 = c1(ios:len(c1))
    ios = index(c1,'"')
    c1 = c1(ios+1:len(c1))
    ios = index(c1,'"')
    c1 = c1(1:ios-1)
    read(c1,*) GPR%cell(i,1,1),GPR%cell(i,2,1),GPR%cell(i,3,1), &
     &     GPR%cell(i,1,2),GPR%cell(i,2,2),GPR%cell(i,3,2), &
     &     GPR%cell(i,1,3),GPR%cell(i,2,3),GPR%cell(i,3,3)
   enddo
  endif
  if (GPR%verbose) write(*,*) 'Got input frame with ',nat,'atoms'

end subroutine

!****************************************************************************************************************
subroutine print_predictions(GPR,un)
 implicit none

  type(SAGPR_Model), intent(inout) :: GPR
  integer j,k,l,un

  if (GPR%atomic) then
   write(un,*) size(GPR%prediction_lm_c,1)
   write(un,*) '# Total',((sum(GPR%prediction_lm_c(:,j,k)),j=1,GPR%degen),k=1,GPR%nw)
  endif
  do l=1,size(GPR%prediction_lm_c,1)
   if (.not. GPR%atomic) then
    write(un,*) ((GPR%prediction_lm_c(l,j,k),j=1,GPR%degen),k=1,GPR%nw)
   else
    write(un,*) GPR%atname_at(l),((GPR%prediction_lm_c(l,j,k),j=1,GPR%degen),k=1,GPR%nw)
   endif
  enddo
  flush(un)

end subroutine

!****************************************************************************************************************

end module
