subroutine phasecomb(nat,nspecies,nneigh,nG,coords,Gvec,phase) 

! This routine computes the Fourier components of the 3D potential by combining the isotropic single-site contributions

implicit none
integer:: iG,i,j,ispe
real*8:: pi
real*8:: coordx,coordy,coordz,arg,Gx,Gy,Gz,cosine,sine,pidiv,pires,sig
integer:: nat,nG,nspecies
integer,dimension(nat,nspecies):: nneigh
real*8,dimension(3,nG):: Gvec
real*8,dimension(nat,nspecies,nat,3):: coords 
complex*16,dimension(2,nspecies,nat,nG):: phase

pi=4.d0*atan(1.d0)
phase(:,:,:,:) = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DEFAULT(private) &
!$OMP SHARED(phase,Gvec,coords,nneigh,nG,nspecies,nat)
!$OMP DO SCHEDULE(dynamic)
do iG=2,nG
   Gx = Gvec(1,iG)
   Gy = Gvec(2,iG)
   Gz = Gvec(3,iG)
   do i=1,nat
      do ispe=1,nspecies
         do j=1,nneigh(i,ispe)
            coordx = coords(i,ispe,j,1)
            coordy = coords(i,ispe,j,2)
            coordz = coords(i,ispe,j,3)
            arg = Gx*coordx + Gy*coordy + Gz*coordz
            cosine = dcos(arg)
            sine = dsin(arg)
            phase(1,ispe,i,iG) = phase(1,ispe,i,iG) + dcmplx(cosine,0.d0)
            phase(2,ispe,i,iG) = phase(2,ispe,i,iG) + dcmplx(0.d0,-sine)
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end
!=================================================================================================================================================

subroutine gcontra(nat,nspecies,nmax,lmax,nG,orthoradint,harmonic,phase,omega) 

! This routine performes the G-contraction needed to get the Poisson-Rayleigh projections <anlm|V_i> 

implicit none
integer:: iG,iat,ispe,n,l,im,lm
integer:: nG,lmax,nmax,nat,nspecies
real*8,dimension(lmax+1,nmax,nG):: orthoradint 
complex*16,dimension(2,nspecies,nat,nG):: phase
complex*16,dimension((lmax+1)*(lmax+1),nG):: harmonic 
complex*16,dimension(nat,nspecies,nmax,lmax+1,2*lmax+1):: omega

!f2py intent(in) nat,nspecies,nmax,lmax,nG,orthoradint,harmonic,phase
!f2py intent(out) omega 
!f2py depend(nat) phase,omega 
!f2py depend(nspecies) phase,omega 
!f2py depend(nmax) orthoradint,omega 
!f2py depend(lmax) orthoradint,harmonic,omega
!f2py depend(nG) orthoradint,harmonic,phase

omega = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DEFAULT(private) &
!$OMP SHARED(nat,nspecies,nmax,lmax,nG,orthoradint,harmonic,phase,omega)
!$OMP DO SCHEDULE(dynamic)
do iat=1,nat
   do iG=2,nG
      do ispe=1,nspecies
         do n=1,nmax
            lm = 1
            do l=1,lmax+1
               if (mod(l-1,2).eq.0) then
                   do im=1,2*(l-1)+1
                      omega(iat,ispe,n,l,im) = omega(iat,ispe,n,l,im) + harmonic(lm,iG) * orthoradint(l,n,iG) * phase(1,ispe,iat,iG)
                      lm = lm + 1
                   end do
                else
                   do im=1,2*(l-1)+1
                      omega(iat,ispe,n,l,im) = omega(iat,ispe,n,l,im) + harmonic(lm,iG) * orthoradint(l,n,iG) * phase(2,ispe,iat,iG)
                      lm = lm + 1
                   end do
                end if 
            end do
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end
