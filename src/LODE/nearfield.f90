subroutine nearfield(nat,nspecies,nmax,lmax,ngrid,nneigh,alpha,coords,sgrid,orthorad,harmonic,weights,omega) 

! This routine performes the G-contraction needed to get the Poisson-Rayleigh projections <anlm|V_i> 

implicit none
integer:: igrid,iat,jat,ispe,n,l,im,lm
integer:: ngrid,lmax,nmax,nat,nspecies
real*8:: alpha,potential,r
integer,dimension(nat,nspecies)::nneigh
real*8,dimension(nat,nspecies,nat,3):: coords
real*8,dimension(nmax,ngrid):: orthorad
real*8,dimension(ngrid):: weights 
real*8,dimension(ngrid,3):: sgrid 
complex*16,dimension((lmax+1)*(lmax+1),ngrid):: harmonic 
complex*16,dimension(nat,nspecies,nmax,lmax+1,2*lmax+1):: omega

omega = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DEFAULT(private) &
!$OMP SHARED(nat,nspecies,nmax,lmax,ngrid,nneigh,alpha,coords,sgrid,orthorad,harmonic,weights,omega)
!$OMP DO SCHEDULE(dynamic)
do iat=1,nat
   do ispe=1,nspecies
      do igrid=1,ngrid
         potential = 0.d0
         do jat=1,nneigh(iat,ispe)
            r = dsqrt((sgrid(igrid,1)-coords(iat,ispe,jat,1))**2 &
                     +(sgrid(igrid,2)-coords(iat,ispe,jat,2))**2 &
                     +(sgrid(igrid,3)-coords(iat,ispe,jat,3))**2)
            potential = potential + erf(dsqrt(alpha)*r)/r
         enddo
         do n=1,nmax
            lm = 1
            do l=1,lmax+1
               do im=1,2*(l-1)+1
                  omega(iat,ispe,n,l,im) = omega(iat,ispe,n,l,im) + harmonic(lm,igrid) &
                                                                  * orthorad(n,igrid)  &
                                                                  * weights(igrid) & 
                                                                  * potential 
                  lm = lm + 1
               end do
            end do
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end

!============================================================================================================================================

subroutine nearfield_ewald(nat,nspecies,nmax,lmax,ngrid,nneigh,neighmax,alpha,coords,sgrid,orthorad,harmonic,weights,sigewald,omega) 

! This routine performes the G-contraction needed to get the Poisson-Rayleigh projections <anlm|V_i> 

implicit none
integer:: igrid,iat,jat,ispe,n,l,im,lm
integer:: ngrid,lmax,nmax,nat,nspecies,neighmax
real*8:: alpha,potential,r,sigewald,alphaewald
integer,dimension(nat,nspecies)::nneigh
!real*8,dimension(3,neighmax,nspecies,nat):: coords
real*8,dimension(nat,nspecies,nat,3):: coords
real*8,dimension(nmax,ngrid):: orthorad
real*8,dimension(ngrid):: weights 
real*8,dimension(ngrid,3):: sgrid 
complex*16,dimension((lmax+1)*(lmax+1),ngrid):: harmonic 
!complex*16,dimension(2*lmax+1,lmax+1,nmax,nspecies,nat):: omega
complex*16,dimension(nat,nspecies,nmax,lmax+1,2*lmax+1):: omega
integer i,j,k

alphaewald = 1.d0/(2.d0*sigewald**2)

	do i=1,nat
	do j=1,nspecies
	write(*,*) 'CHECK NNEIGH',nneigh(i,j)
	enddo
	enddo
	do i=1,nat
	do j=1,nspecies
	do k=1,nat
	do l=1,3
	write(*,*) 'CHECK COORDS',coords(i,j,k,l)
	enddo
	enddo
	enddo
	enddo
	do i=1,nmax
	do j=1,ngrid
	write(*,*) 'CHECK ORTHORAD',orthorad(i,j)
	enddo
	enddo
	do i=1,ngrid
	write(*,*) 'CHECK WEIGHTS',weights(i)
	enddo
	do i=1,ngrid
	do j=1,3
	write(*,*) 'CHECK SGRID',sgrid(i,j)
	enddo
	enddo
	do i=1,(lmax+1)*(lmax+1)
	do j=1,ngrid
	write(*,*) 'CHECK HARMONIC',harmonic(i,j)
	enddo
	enddo
	write(*,*) nat,nspecies,neighmax

omega = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DEFAULT(private) &
!$OMP SHARED(nat,nspecies,nmax,lmax,ngrid,nneigh,alpha,alphaewald,coords,sgrid,orthorad,harmonic,weights,omega)
!$OMP DO SCHEDULE(dynamic)
do iat=1,nat
   do ispe=1,nspecies
      do igrid=1,ngrid
         potential = 0.d0
         do jat=1,nneigh(iat,ispe)
!            r = dsqrt((sgrid(igrid,1)-coords(1,jat,ispe,iat))**2 &
!                     +(sgrid(igrid,2)-coords(2,jat,ispe,iat))**2 &
!                     +(sgrid(igrid,3)-coords(3,jat,ispe,iat))**2)
            r = dsqrt((sgrid(igrid,1)-coords(iat,ispe,jat,1))**2 &
                     +(sgrid(igrid,2)-coords(iat,ispe,jat,2))**2 &
                     +(sgrid(igrid,3)-coords(iat,ispe,jat,3))**2)
            potential = potential + ( erf(dsqrt(alpha)*r) - erf(dsqrt(alphaewald)*r) ) / r
         enddo
         do n=1,nmax
            lm = 1
            do l=1,lmax+1
               do im=1,2*(l-1)+1
                  omega(iat,ispe,n,l,im) = omega(iat,ispe,n,l,im) + harmonic(lm,igrid) &
!                  omega(im,l,n,ispe,iat) = omega(im,l,n,ispe,iat) + harmonic(lm,igrid) &
                                                                  * orthorad(n,igrid)  &
                                                                  * weights(igrid) & 
                                                                  * potential 
                  lm = lm + 1
               end do
            end do
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end

