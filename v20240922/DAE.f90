!==================================================================================================
module DAE

   use, intrinsic :: iso_c_binding
   use fsundials_core_mod

contains

! ----------------------------------------------------------------
! DAE residual function
! ----------------------------------------------------------------
integer(c_int) function residuals(time, sunvec_y, sunvec_yp, sunvec_r, user_data) result(ierr) bind(C,name='residuals')

   use fnvector_serial_mod
   use globals
 
   implicit none
   real(c_double), value :: time   ! current time
   type(N_Vector) :: sunvec_y      ! solution N_Vector
   type(N_Vector) :: sunvec_yp     ! derivative N_Vector
   type(N_Vector) :: sunvec_r      ! residual N_Vector
   type(c_ptr), value :: user_data ! user-defined data
 
!  pointers to data in SUNDIALS vectors
   real(c_double), pointer :: y(:)
   real(c_double), pointer :: y_t(:)
   real(c_double), pointer :: r(:)
 
   integer i,j,n
   real(c_double) polp,fgr

   real(c_double) rr(maxeq)

   real(c_double) tmp(maxtab)
  
!  get data arrays from SUNDIALS vectors
   y => FN_VGetArrayPointer(sunvec_y)
   y_t => FN_VGetArrayPointer(sunvec_yp)
   r => FN_VGetArrayPointer(sunvec_r)
 
   call read_from_y(y)
   call read_from_y_t(y_t)

   do j=1,nzz
!     given deformations, update geometry
      call update_geom(j)
!     linear power (W/m)
      ql(j) = qqv1(j)*azf0
   end do

!  counter of equations
   n = 0
!  residuals for whole fuel rod   
   call residuals_0(time, rr)
   do i = 1, 3
      n = n + 1
      r(n) = rr(i)
   end do
   do j=1,nzz
!     residuals for axial layer j
      call residuals_j(j, time, rr)
      do i = 1, neq_j
         n = n + 1
         r(n) = rr(i)
      end do
   end do
 
!  return success
   ierr = 0
   return
   
end function residuals

! ----------------------------------------------------------------
! Returns residual r for the whole fuel rod
! ----------------------------------------------------------------
subroutine residuals_0(time, r)

   use globals
   implicit none
   integer j, n
   real(c_double), intent(inout) :: r(:)
   real(c_double) time,fgGrate,fgRrate,vol,rate
   real(c_double) fgr,gpresf

   n = 0
!  fission gas generation rate (mol/s)
   fgGrate = 0.0d0
!  fission gas release rate (mol/s)
   fgRrate = 0.0d0
   do j=1,nzz
!     volume of fuel axial layer
      vol = azf0*dz0(j)
!     fission rate (fiss/s): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
      rate = qqv1(j)*vol/(1.60214d-13*210.0d0)
!     fission gas generation rate (mol/s) 25 atoms of FGP are produced per 100 fissions (Olander)
      rate = 0.25d0*rate / avo
      fgGrate = fgGrate + rate
!     fission gas release rate (mol/s)
      fgRrate = fgRrate + rate * fgr(j)
   end do
   if(fgGrate .lt. 0.0d0)then
      fgGrate = 0.0d0
      fgRrate = 0.0d0
   end if

!  fission gas production rate ODE (mol/s)
   n = n + 1
   r(n) = dfggen - fgGrate
   
!  fission gas release rate ODE (mol/s)
   n = n + 1
   r(n) = dfgrel - fgRrate
 
!  gas pressure AE (MPa)
   n = n + 1
   r(n) = gpres - gpresf(time)

end
! ----------------------------------------------------------------
! Returns residual r for axial layer j
! ----------------------------------------------------------------
subroutine residuals_j(j, time, r)

   use globals
   implicit none
   real(c_double) time
   real(c_double), intent(inout) :: r(:)

   integer j, i, n
   real(c_double) polp,clamb,ctexp,gaphtc,flamb,ftexp,felmod,fpoir,fcp,ccp,celmod,cpoir,fgr,gpresf,fswel, &
  &               fcreep,ccreep
   real(c_double) tmp(maxtab),kclad_,kfuel_,sto,cp,vol, &
  &       tem_(maxr),qf_(maxr), &
  &       fyng(maxr),fpoi(maxr), &
  &       efr_(maxr),sigfh_(maxr), &
  &       cyng(maxr),cpoi(maxr), &
  &       er_(maxr),sigh_(maxr), &
  &       eel

   n = 0
   do i=1,nf
!     fuel young modulus and poisson ratio
      fyng(i) = felmod(tem(i,j),rof0,pucont,fmat)
      fpoi(i) = fpoir(fmat)
!     effective fuel stress
      sigf(i,j) = dsqrt((sigfh(i,j)-sigfz(i,j))**2 + &
     &                  (sigfh(i,j)-sigfr(i,j))**2 + &
     &                  (sigfz(i,j)-sigfr(i,j))**2) / dsqrt(2.0d0)       
   end do
   do i=1,nc
!     clad young modulus and poisson ratio
      cyng(i) = celmod(tem(nf+i,j),cmat)
      cpoi(i) = cpoir(cmat)
!     effective clad stress
      sig(i,j) = dsqrt((sigh(i,j)-sigz(i,j))**2 + &
     &                 (sigh(i,j)-sigr(i,j))**2 + &
     &                 (sigz(i,j)-sigr(i,j))**2) / dsqrt(2.0d0)
   end do
   do i=1,nf-1
      efr_(i) = 0.5d0*(efr(i,j) + efr(i+1,j))
      sigfh_(i) = 0.5d0*(sigfh(i,j) + sigfh(i+1,j))
   end do
   do i=1,nc-1
      er_(i) = 0.5d0*(er(i,j) + er(i+1,j))
      sigh_(i) = 0.5d0*(sigh(i,j) + sigh(i+1,j))
   end do

!  fuel power density AE (W/m3)
   do i=1,nqv
      tmp(i) = qv(j,i)
   end do
   n = n + 1
   r(n) = qqv1(j) - polp(time,tmp,tqv,nqv)

!  fuel burnup rate ODE (J/kg-s)
   n = n + 1
   r(n) = dbup(j) - qqv1(j)/rof0
 
!  gap width AE
   n = n + 1
   r(n) = gap(j) - (rci(j) - rfo(j))
   
!  gap state (flag)
   if(gap(j) .le. (ruff + rufc))flag(j) = 'clos'

!   the condition of the gap reopening was commented -- to be further explored
!   if(flag(j) .eq. 'clos' .and. pfc(j) .lt. gpres .and. dpfc(j) .lt. 0.0d0)flag(j) = 'open'

!  contact pressure AE (MPa)
   n = n + 1
   if(flag(j) .eq. 'clos')then
      r(n) = pfc(j) - dabs(sigr(1,j))
   else      
      r(n) = pfc(j)
   end if

!  gap conductance AE (W/m2K)
   n = n + 1
   r(n) = hgapt(j) - gaphtc(j)
   
!  calculate temperatures at boundaries between nodes
   do i=1,nf+nc-1
      tem_(i) = 0.5d0*(tem(i,j) + tem(i+1,j))
   end do

!  table of fuel stoichiometry vs fuel burnup (MWd/kg)
   if(nfsto.gt.0)then
      sto = polp(bup(j),stob,bsto,nfsto)
   else
      sto = sto0
   end if
   
!  heat fluxes inside fuel pellet
   do i=1,nf-1
      kfuel_ = flamb(rof0,bup(j),tem_(i),pucont,sto,fmat)
      qf_(i) = kfuel_*(tem(i,j)-tem(i+1,j))/drf0
   end do

!  heat flux between fuel pellet and cladding
   qf_(nf) = hgapt(j)*(tem(nf,j) - tem(nf+1,j))
!  heat fluxes inside cladding
   do i=nf+1,nf+nc-1
      kclad_=clamb(tem_(i),cmat)
      qf_(i) = kclad_*(tem(i,j)-tem(i+1,j))/drc0
   end do

!  fuel temperature ODE
   do i=1,nf
      cp = fcp(tem(i,j),pucont,fmat)
      n = n + 1
      r(n) = qqv1(j)*az0(i) - qf_(i)*2.0d0*pi*rad_0(i) - rof0*cp*az0(i)*dtem(i,j)
      if(i .gt. 1) r(n) = r(n) + qf_(i-1)*2.0d0*pi*rad_0(i-1)
   end do

!  clad temperature ODE
   do i=nf+1,nf+nc-1
      cp = ccp(tem(i,j),cmat)
      n = n + 1
      r(n) = qf_(i-1)*2.0d0*pi*rad_0(i-1) - qf_(i)*2.0d0*pi*rad_0(i) - roc0*cp*az0(i)*dtem(i,j)
   end do

!  outer cladding temperature AE
   do i=1,ntco
      tmp(i)=tco(j,i)
   end do
   n = n + 1
   r(n) = tem(nf+nc,j) - polp(time,tmp,ttco,ntco)

   if(inomech .eq. 0) then

!     FUEL: swelling rate ODE
      do i=1,nf
         n = n + 1
         r(n) = defs(i,j) - fswel(i,j)
      end do
      
!     FUEL: thermal expansion AE
      do i=1,nf
         n = n + 1
         r(n) = eft(i,j) - (ftexp(tem(i,j),pucont,fmat) - ftexp(293.15d0,pucont,fmat))
      end do
      
!     FUEL: creep rate ODEs
      if(ifcreep .eq. 1)then
!        fuel creep effective strain rate ODE
         do i=1,nf
            n = n + 1
            r(n) = defce(i,j) - fcreep(i,j)
         end do
      
!        fuel creep hoop strain rate ODE: Prandtl-Reuss flow rule
         do i=1,nf
            n = n + 1
            if(sigf(i,j) .eq. 0.0d0)then
               r(n) = defch(i,j)
            else
               r(n) = defch(i,j) - 1.5d0*defce(i,j)*( sigfh(i,j) - (sigfh(i,j) + sigfr(i,j) + sigfz(i,j))/3.0d0 )/sigf(i,j)
            end if
         end do
         
!        FUEL: creep radial strain rate ODE: Prandtl-Reuss flow rule
         do i=1,nf
            n = n + 1
            if(sigf(i,j) .eq. 0.0d0)then
               r(n) = defcr(i,j)
            else
               r(n) = defcr(i,j) - 1.5d0*defce(i,j)*( sigfr(i,j) - (sigfh(i,j) + sigfr(i,j) + sigfz(i,j))/3.0d0 )/sigf(i,j)
            end if
         end do
         
!        FUEL: creep axial strain rate ODE: Prandtl-Reuss flow rule
         do i=1,nf
            n = n + 1
            if(sigf(i,j) .eq. 0.0d0)then
               r(n) = defcz(i,j)
            else
               r(n) = defcz(i,j) - 1.5d0*defce(i,j)*( sigfz(i,j) - (sigfh(i,j) + sigfr(i,j) + sigfz(i,j))/3.0d0 )/sigf(i,j)
            end if
         end do
      end if
      
!     FUEL: Hooke's eq. for eh:yng*eh-sigh+pss*sigr+pss*sigz=yng*eh_pct
      do i=1,nf
         n = n + 1
         eel = efh(i,j) - eft(i,j) - efs(i,j)/3.0d0 - efch(i,j) !- efd(i,j)
         r(n) = fyng(i)*eel - sigfh(i,j) + fpoi(i)*(sigfr(i,j) + sigfz(i,j))
      end do
      
!     FUEL: Hooke's eq. for er:yng*er+pss*sigh-sigr+pss*sigz=yng*er_pct
      do i=1,nf
         n = n + 1
         eel = efr(i,j) - eft(i,j) - efs(i,j)/3.0d0 - efcr(i,j) !- efd(i,j)
         r(n) = fyng(i)*eel - sigfr(i,j) + fpoi(i)*(sigfh(i,j) + sigfz(i,j))
      end do
      
!     FUEL: Hooke's eq. for ez:yng*ez+pss*sigh+pss*sigr-sigz=yng*ez_pct
      do i=1,nf
         n = n + 1
         eel = efz(j) - eft(i,j) - efs(i,j)/3.0d0 - efcz(i,j) !- efd(i,j)
         r(n) = fyng(i)*eel - sigfz(i,j) + fpoi(i)*(sigfh(i,j) + sigfr(i,j))
      end do
      
!     FUEL: strain compatibility eq.:(rr-rl)*er+rl0*ehl-rr0*ehr=0
      do i=1,nf-1
         n = n + 1
         r(n) = drf(i,j)*efr_(i) + rad(i,j)*efh(i,j) - rad(i+1,j)*efh(i+1,j)
      end do
      
!     FUEL: stress equilibrity eq.:(rr-rl)*sigh+rl0*sigrl-rr0*sigrr=0
      do i=1,nf-1
         n = n + 1
         r(n) = drf(i,j)*sigfh_(i) + rad(i,j)*sigfr(i,j) - rad(i+1,j)*sigfr(i+1,j)
      end do
      
!     FUEL: boundary conditions - 1
      n = n + 1
      r(n) = 0.0d0
      if(flag(j).eq.'clos')then
         do i=1,nf
            r(n) = r(n) + az0(i)*sigfz(i,j)
         end do
         do i=1,nc
            r(n) = r(n) + az0(nf+i)*sigz(i,j)
         end do
         r(n) = r(n) + pcool*pi*rco0**2
      else
         do i=1,nf
            r(n) = r(n) + az0(i)*sigfz(i,j)
         end do
         r(n) = r(n) + azf0*gpres
      end if
      
!     FUEL: boundary conditions - 2
      n = n + 1
      if(rad0(1).eq.0.d0)then
!        no central hole. symmetry: sigr=sigh
         r(n) = sigfr(1,j) - sigfh(1,j)
      else
!        central hole. gas pressure: sigr=-gpres
         r(n) = sigfr(1,j) + gpres
      end if
      
!     FUEL: boundary conditions - 3
      n = n + 1
      if(flag(j).eq.'clos')then
!        closed gap: D(ez_fuel)=D(ez_clad)
         r(n) = defz(j) - dez(j)
      else
!        open gap: sigr=-gpres
         r(n) = sigfr(nf,j) + gpres
      end if
      
!     CLAD: thermal expansion rate AE
      do i=1,nc
         n = n + 1
         r(n) = et(i,j) - (ctexp(tem(nf+i,j),cmat) - ctexp(293.15d0,cmat))
      end do
      
!     CLAD: creep rate ODEs
      if(iccreep .eq. 1)then
!        clad creep effective strain rate ODE (1/s)
         do i=1,nc
            n = n + 1
            r(n) = dece(i,j) - ccreep(tem(nf+i,j),sig(i,j),cmat)
         end do
      
!        clad creep hoop strain rate ODE (1/s): Prandtl-Reuss flow rule
         do i=1,nc
            n = n + 1
            if(sig(i,j) .eq. 0.0d0)then
               r(n) = dech(i,j)
            else
               r(n) = dech(i,j) - 1.5d0*dece(i,j)*( sigh(i,j) - (sigh(i,j) + sigr(i,j) + sigz(i,j))/3.0d0 )/sig(i,j)
            end if
         end do
         
!        clad creep radial strain rate ODE (1/s): Prandtl-Reuss flow rule
         do i=1,nc
            n = n + 1
            if(sig(i,j) .eq. 0.0d0)then
               r(n) = decr(i,j)
            else
               r(n) = decr(i,j) - 1.5d0*dece(i,j)*( sigr(i,j) - (sigh(i,j) + sigr(i,j) + sigz(i,j))/3.0d0 )/sig(i,j)
            end if
         end do
         
!        clad creep axial strain rate ODE (1/s): Prandtl-Reuss flow rule
         do i=1,nc
            n = n + 1
            if(sig(i,j) .eq. 0.0d0)then
               r(n) = decz(i,j)
            else
               r(n) = decz(i,j) - 1.5d0*dece(i,j)*( sigz(i,j) - (sigh(i,j) + sigr(i,j) + sigz(i,j))/3.0d0 )/sig(i,j)
            end if
         end do
      end if
      
!     CLAD: Hooke's eq. for eh:yng*eh-sigh+pss*sigr+pss*sigz=yng*eh_pct
      do i=1,nc
         n = n + 1
         eel = eh(i,j)-et(i,j)-ech(i,j)
         r(n) = cyng(i)*eel - sigh(i,j) + cpoi(i)*(sigr(i,j) + sigz(i,j))
      end do
      
!     CLAD: Hooke's eq. for er:yng*er+pss*sigh-sigr+pss*sigz=yng*er_pct
      do i=1,nc
         n = n + 1
         eel = er(i,j)-et(i,j)-ecr(i,j)
         r(n) = cyng(i)*eel - sigr(i,j) + cpoi(i)*(sigh(i,j) + sigz(i,j))
      end do
      
!     CLAD: Hooke's eq. for ez:yng*ez+pss*sigh+pss*sigr-sigz=yng*ez_pct
      do i=1,nc
         n = n + 1
         eel = ez(j)-et(i,j)-ecz(i,j)
         r(n) = cyng(i)*eel - sigz(i,j) + cpoi(i)*(sigh(i,j) + sigr(i,j))
      end do
      
!     CLAD: strain compatibility eq.:(rr-rl)*er+rl0*ehl-rr0*ehr=0
      do i=1,nc-1
         n = n + 1
         r(n) = drc(i,j)*er_(i) + rad(nf+i,j)*eh(i,j) - rad(nf+i+1,j)*eh(i+1,j)
      end do
      
!     CLAD: stress equilibrity eq.:(rr-rl)*sigh+rl0*sigrl-rr0*sigrr=0
      do i=1,nc-1
         n = n + 1
         r(n) = drc(i,j)*sigh_(i) + rad(nf+i,j)*sigr(i,j) - rad(nf+i+1,j)*sigr(i+1,j)
      end do
      
!     CLAD: boundary conditions - 1
      n = n + 1
      if(flag(j).eq.'clos')then
!        closed gap: sigr_fuel=sigr_clad
         r(n) = sigfr(nf,j) - sigr(1,j)
      else
!        open gap: axial stress equilibrity eq.
         r(n) = 0.0d0
         do i=1,nc
            r(n) = r(n) + az0(nf+i)*sigz(i,j)
         end do
         r(n) = r(n) + pcool*pi*rco0**2 - gpres*pi*rci0**2
      end if
      
!     CLAD: boundary conditions - 2
      n = n + 1
      if(flag(j).eq.'clos')then
!        closed gap: D(eh_fuel)=D(eh_clad)
         !r(n) = defh(nf,j) - deh(1,j)
         r(n) = gap(j) - (ruff+rufc)
      else
!        open gap: gas pressure: sigr=-gpres
         r(n) = sigr(1,j) + gpres
      end if
      
!     CLAD: boundary conditions - 3
      n = n + 1
!     sigr=-pcool      
      r(n) = sigr(nc,j) + pcool
   end if

end

! ----------------------------------------------------------------
! Find index of the given equation id in the structure id_eq
! ----------------------------------------------------------------
integer function indx(str,j,i)
   use globals
   implicit none
   integer j,i,n
   character (len=10) str

   indx = 0
   do n = 1,neq
      if(id_eq(n)%id .eq. str .and. &
     &   id_eq(n)%j .eq. j .and. &
     &   id_eq(n)%i .eq. i) indx = n
   end do
end function

! ----------------------------------------------------------------
! Root finding function
! ----------------------------------------------------------------
integer(c_int) function root(t, sunvec_y, sunvec_yp, gout, user_data) result(ierr) bind(C,name='root')

   use, intrinsic :: iso_c_binding
   use fnvector_serial_mod
   use globals

   implicit none

   real(c_double), value :: t       ! current time
   type(N_Vector) :: sunvec_y       ! solution N_Vector
   type(N_Vector) :: sunvec_yp      ! derivative N_Vector
   real(c_double) :: gout(nzz)      ! root function values
   type(c_ptr), value :: user_data  ! user-defined data

!  pointers to data in SUNDIALS vectors
   real(c_double), pointer :: y(:)

   integer j,k
    
!  get data arrays from SUNDIALS vectors
   y => FN_VGetArrayPointer(sunvec_y)

   call read_from_y(y)

!  fill root vector
   k = 0
   do j = 1,nzz
      k = k + 1
      if(gap(j) .gt. (ruff + rufc) .and. flag(j) .eq. 'open')then
         gout(k) = dmax1(0.0d0, gap(j) - (ruff + rufc) )
      else
         gout(k) = -1
      end if
   end do

!  return success
   ierr = 0
   return

end function root

end module DAE
