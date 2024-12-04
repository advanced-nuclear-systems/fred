!--------------------------------------------------------------------------------------------------
module FRED

   use, intrinsic :: iso_c_binding
   use fsundials_core_mod
   use globals
   implicit none

contains

!--------------------------------------------------------------------------------------------------
! DAE residual function
!--------------------------------------------------------------------------------------------------
   integer(c_int) function residuals(time, sunvec_y, sunvec_yp, sunvec_r, user_data) result(ierr) bind(C,name='residuals')

   use fnvector_serial_mod
 
   real(c_double), value :: time   ! current time
   type(N_Vector) :: sunvec_y      ! solution N_Vector
   type(N_Vector) :: sunvec_yp     ! derivative N_Vector
   type(N_Vector) :: sunvec_r      ! residual N_Vector
   type(c_ptr), value :: user_data ! user-defined data
 
!  pointers to data in SUNDIALS vectors
   real(c_double), pointer :: y(:)
   real(c_double), pointer :: y_t(:)
   real(c_double), pointer :: r(:)
 
   integer(c_int) i,j,n
   real(c_double) fgr

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

!--------------------------------------------------------------------------------------------------
! Returns residual r for the whole fuel rod
!--------------------------------------------------------------------------------------------------
   subroutine residuals_0(time, r)

   integer(c_int) j, n
   real(c_double), intent(inout) :: r(:)
   real(c_double) time,fgGrate,fgRrate,vol,rate
   real(c_double) gpresf

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
!--------------------------------------------------------------------------------------------------
! Returns residual r for axial layer j
!--------------------------------------------------------------------------------------------------
   subroutine residuals_j(j, time, r)

   real(c_double) time
   real(c_double), intent(inout) :: r(:)

   integer(c_int) j, i, n
   real(c_double) polp,polp2,clamb,ctexp,gaphtc,flamb_,ftexp,felmod_,fpoir,fcp,ccp,celmod,cpoir,gpresf,fswel, &
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
      fyng(i) = felmod(tem(i,j),rof0,pucont)
      fpoi(i) = fpoir()
!     effective fuel stress
      sigf(i,j) = dsqrt((sigfh(i,j)-sigfz(i,j))**2 + &
     &                  (sigfh(i,j)-sigfr(i,j))**2 + &
     &                  (sigfz(i,j)-sigfr(i,j))**2) / dsqrt(2.0d0)       
   end do
   do i=1,nc
!     clad young modulus and poisson ratio
      cyng(i) = celmod(tem(nf+i,j))
      cpoi(i) = cpoir()
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
   !r(n) = qqv1(j) - polp2(time,tqv,tmp,nqv,360.0d0)

!  fuel burnup rate ODE (J/kg-s)
   n = n + 1
   r(n) = dbup(j) - qqv1(j)/rof0
 
!  gap width AE
   n = n + 1
   r(n) = gap(j) - (rci(j) - rfo(j))
   
!  gap state (flag)
   if(gap(j) .le. (ruff + rufc)) flag(j) = 'clos'

!   the condition of the gap reopening was commented -- to be further explored
!    if(sigr(1,j) .ge. 0.d0) flag(j) = 'open'

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
      kfuel_ = flamb(rof0,bup(j),tem_(i),pucont,sto)
      qf_(i) = kfuel_*(tem(i,j)-tem(i+1,j))/drf0
   end do

!  heat flux between fuel pellet and cladding
   qf_(nf) = hgapt(j)*(tem(nf,j) - tem(nf+1,j))
!  heat fluxes inside cladding
   do i=nf+1,nf+nc-1
      kclad_=clamb(tem_(i))
      qf_(i) = kclad_*(tem(i,j)-tem(i+1,j))/drc0
   end do

!  fuel temperature ODE
   do i=1,nf
      cp = fcp(tem(i,j),pucont)
      n = n + 1
      r(n) = qqv1(j)*az0(i) - qf_(i)*2.0d0*pi*rad_0(i) - rof0*cp*az0(i)*dtem(i,j)
      if(i .gt. 1) r(n) = r(n) + qf_(i-1)*2.0d0*pi*rad_0(i-1)
   end do

!  clad temperature ODE
   do i=nf+1,nf+nc-1
      cp = ccp(tem(i,j))
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
         r(n) = eft(i,j) - (ftexp(tem(i,j),pucont) - ftexp(293.15d0,pucont))
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
         r(n) = et(i,j) - (ctexp(tem(nf+i,j)) - ctexp(293.15d0))
      end do
      
!     CLAD: creep rate ODEs
      if(iccreep .eq. 1)then
!        clad creep effective strain rate ODE (1/s)
         do i=1,nc
            n = n + 1
            r(n) = dece(i,j) - ccreep(tem(nf+i,j),sig(i,j))
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

!--------------------------------------------------------------------------------------------------
! Find index of the given equation id in the structure id_eq
!--------------------------------------------------------------------------------------------------
   integer(c_int) function indx(str,j,i)

   integer(c_int) j,i,n
   character (len=10) str

   indx = 0
   do n = 1,neq
      if(id_eq(n)%id .eq. str .and. &
     &   id_eq(n)%j .eq. j .and. &
     &   id_eq(n)%i .eq. i) indx = n
   end do
   end function indx

!--------------------------------------------------------------------------------------------------
! Root finding function
!--------------------------------------------------------------------------------------------------
   integer(c_int) function root(t, sunvec_y, sunvec_yp, gout, user_data) result(ierr) bind(C,name='root')

   use fnvector_serial_mod

   real(c_double), value :: t       ! current time
   type(N_Vector) :: sunvec_y       ! solution N_Vector
   type(N_Vector) :: sunvec_yp      ! derivative N_Vector
   real(c_double) :: gout(nzz)      ! root function values
   type(c_ptr), value :: user_data  ! user-defined data

!  pointers to data in SUNDIALS vectors
   real(c_double), pointer :: y(:)

   integer(c_int) j,k
    
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

   do j = 1,ntfix
      k = k + 1
      gout(k) = t - tfix(j)
   end do

!  return success
   ierr = 0
   return

   end function root

!--------------------------------------------------------------------------------------------------
! Calculation of the clad specific heat J/(kg*K)
!     input arguments :
!        tk .... temperature (k)
!--------------------------------------------------------------------------------------------------
   real(c_double) function ccp(tk) bind(C,name='ccp_')

   real(c_double) tk,temp(16),cp(16),temp2(17),cp2(17),tf

   if(cmat.eq.'aim1')then
!     L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
      ccp = 431.0d0 + 0.177d0*tk + 8.72d-5*tk**2
   else if(cmat.eq.'t91')then
      ccp = 431.0d0 + 0.177d0*tk + 8.72d-5*tk**2
   else
      write(*,*)'ccp: wrong clad material:',cmat
      stop
   end if
   end function ccp

!--------------------------------------------------------------------------------------------------
! Returns clad material effective creep rate (1/s)
!     input arguments :
!        tk -- temperature (K)
!        sg -- effective stress (MPa)
!--------------------------------------------------------------------------------------------------
   real(c_double) function ccreep(tk,sg) bind(C,name='ccreep_')

   real(c_double) tk,sg,tt,E,nflux

   if(cmat.eq.'aim1')then
!     L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
      tt = 1.986d0*tk
!     thermal creep (%/h)
      ccreep = 2.3d14*dexp(-8.46d4/tt)*dsinh(39.72d0*sg/tt)
!     mean neutron energy (MeV)
      E = 2.0d0
!     neutron flux (n/cm2-s)
      nflux = 1.0d18
!     irradiation creep (%/h)
      ccreep = ccreep + 3.2d-24*E*nflux*sg
!     convert from %/h to 1/s
      ccreep=ccreep/100.d0/3.6d3
   else
      ccreep=0.0d0
   end if
   return
   end function ccreep

!--------------------------------------------------------------------------------------------------
! Calculation of the clad Young's modulus
!                unit: MPa
!     input arguments :
!        tk .... temperature (K)
!--------------------------------------------------------------------------------------------------
   real(c_double) function celmod(tk) bind(C,name='celmod_')

   real(c_double) tk,tc

   celmod = 0.0d0
   if(cmat.eq.'aim1')then
      tc=tk-273.15d0
      celmod=2.027d11-8.167d7*tc
   else if(cmat.eq.'t91')then
      tc=tk-273.15
      if(tc .le. 500.0)then
         celmod=2.073d11-64.58d6*tc
      else if(tc .le. 600.0)then
         celmod=2.95d11-240.0d6*tc
      end if   
   else
      write(*,*)'celmod: wrong clad material:',cmat
      stop
   end if
   celmod = celmod/1.0d6
   return
   end function celmod

!--------------------------------------------------------------------------------------------------
! Calculation of the clad thermal conductivity
!                unit: w/(m*k)
!     input arguments :
!        tk .... temperature (K)
!--------------------------------------------------------------------------------------------------
   real(c_double) function clamb(tk) bind(C,name='clamb_')

   real(c_double) tk, tc

   if(cmat.eq.'aim1')then
      tc=tk-273.15d0
      clamb=13.95d0+1.163d-2*tc
   else if(cmat.eq.'t91')then
      tc=tk-273.15d0
      clamb = 23.71 + 0.01718*tc - 1.45E-05*tc**2
   else
      write(*,*)'clamb: wrong clad material:',cmat
      stop
   end if
   return
   end function clamb

!--------------------------------------------------------------------------------------------------
! Calculation of the clad Poisson ratio
!                unitless
!     input arguments :
!        tk .... temperature (k)
!--------------------------------------------------------------------------------------------------
   real(c_double) function cpoir() bind(C,name='cpoir_')

   if(cmat.eq.'aim1')then
      cpoir=0.289d0
   else if(cmat.eq.'t91')then
      cpoir=0.28d0
   else
      write(*,*)'cpoir: wrong clad material:',cmat
      stop
   end if
   return
   end function cpoir

!--------------------------------------------------------------------------------------------------
! Calculates burst stress for cladding
!--------------------------------------------------------------------------------------------------
   real(c_double) function csigb(tk) bind(C,name='csigb_')
   real(c_double) tk

   if(cmat.eq.'aim1'.or.cmat.eq.'t91')then
      csigb = 1.5957d9 - 4.7253d6*tk + 9.8851d3*tk**2 - 8.8864d0*tk**3 + 2.5538d-3*tk**4
   else
      csigb = 0.0
   end if
   return
   end function csigb

!--------------------------------------------------------------------------------------------------
! Calculate yield stress for cladding
!--------------------------------------------------------------------------------------------------
   real(c_double) function csigy(tk) bind(C,name='csigy_')
   real(c_double) tk, tc

   if(cmat.eq.'t91')then
      csigy = 1.3109d9 - 3.6916d6*tk + 7.8909d3*tk**2 - 7.3551d0*tk**3 + 2.1966d-3*tk**4
   else if(cmat.eq.'aim1')then
     tc = tk - 273.15
     if(tc.lt.600)then
         csigy = 5.555d8 - 0.25d0*tc
     else if(tc.ge.600 .and. tc.le.1000)then
         csigy = 4.055d8 - 0.755d0*(tc-600.0d0)
     else
         csigy = 3.455d8 - 0.25d0*tc
     end if
   else
      csigy = 0.0
   end if
   return
   end function csigy

!--------------------------------------------------------------------------------------------------
! Calculates the clad thermal expansion strain
!--------------------------------------------------------------------------------------------------
   real(c_double) function ctexp(tk) bind(C,name='ctexp_')

   real(c_double) tk, tc

   if(cmat.eq.'aim1')then
!     L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
      tc = tk - 273.15d0
      ctexp=-3.101d-4 + 1.545d-5*tc + 2.75d-9*tc**2
      ctexp=-0.2177d-2 + 6.735d-6*tk + 5.12d-9*tk**2 - 2.248d-12*tk**3 + 3.933d-16*tk**4
   else if(cmat.eq.'t91')then
      ctexp=-0.2177d-2 + 6.735d-6*tk + 5.12d-9*tk**2 - 2.248d-12*tk**3+3.933d-16*tk**4
   else
      write(*,*)'wrong cmat option ',cmat
      stop
   end if
   return
   end function ctexp

!--------------------------------------------------------------------------------------------------
! Calculates ultimate elongation for cladding
!--------------------------------------------------------------------------------------------------
   real(c_double) function cuelon(tk) bind(C,name='cuelon_')

   real*8 tk

   if(cmat.eq.'aim1'.or.cmat.eq.'t91')then
      if(tk.le.720.d0)then
         cuelon=-0.58258 + 8.4018d-3*tk - 3.2807d-5*tk**2 + 4.9989d-8*tk**3 - 2.6347d-11*tk**4
      else
         cuelon=-0.85401 + 4.5753d-3*tk - 8.2202d-6*tk**2 + 6.1983d-9*tk**3 - 1.6897d-12*tk**4
      end if           
   else
      cuelon = 0.0
   end if
   return
   end function cuelon

!--------------------------------------------------------------------------------------------------
! Calculation of fuel specific heat
!                  unit : j/(kg*k)
!       input arguments :
!          tk .... temperature (k)
!--------------------------------------------------------------------------------------------------
   real(c_double) function fcp(tk,pucont) bind(C,name='fcp_')

   real(c_double) polp
   real(c_double) tk,pucont,ucp_,pucp_,tPopov(29),cpuo2Popov(29),cppuo2Popov(29),tliq,tf

!  S.G. Popov, J.J. Carbajo, V.K. Ivanov, G.L. Yoder. Thermophysical properties of MOX and UO2 fuels including the effects of irradiation. ORNL/TM-2000/351 http://web.ornl.gov/~webworks/cpr/v823/rpt/109264.pdf
   data tPopov     /300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0, &
  &                1000.d0,1100.d0,1200.d0,1300.d0,1400.d0,1500.d0,1600.d0, &
  &                1700.d0,1800.d0,1900.d0,2000.d0,2100.d0,2200.d0,2300.d0, &
  &                2400.d0,2500.d0,2600.d0,2700.d0,2800.d0,2900.d0,3000.d0, &
  &                3100.d0/
   data cpuo2Popov /235.51d0,265.79d0,282.14d0,292.21d0,299.11d0,304.24d0,308.32d0, &
  &                 311.74d0,314.76d0,317.59d0,320.44d0,323.60d0,327.41d0,332.31d0, &
  &                 338.77d0,347.30d0,358.40d0,372.54d0,390.12d0,411.46d0,436.78d0, &
  &                 466.21d0,499.80d0,537.50d0,579.17d0,624.63d0,673.63d0,725.88d0, &
  &                 781.08d0/
   data cppuo2Popov/203.71d0,225.79d0,235.77d0,241.33d0,244.92d0,247.47d0,249.44d0, &
  &                 251.04d0,252.43d0,253.70d0,254.96d0,256.32d0,257.93d0,259.93d0, &
  &                 262.47d0,265.67d0,269.57d0,274.17d0,279.35d0,284.96d0,290.80d0, &
  &                 296.70d0,302.47d0,308.00d0,313.21d0,318.08d0,322.59d0,326.76d0, &
  &                 330.64d0/

   if(fmat.eq.'mox')then
      ucp_=polp(tk,cpuo2Popov,tPopov,29)
      pucp_=polp(tk,cppuo2Popov,tPopov,29)
      fcp=ucp_*(1.d0-pucont)+pucp_*pucont
!     melting point
      tliq=3120.0d0-388.1d0*pucont-30.4d0*pucont**2
      if(tk.ge.tliq-1.0d0)then
         fcp=fcp+1.d5*abs(tk-tliq+1.0d0)
      end if
   else
      write(*,*)'wrong fuel material:',fmat
      stop
   end if
   return
   end function fcp

!--------------------------------------------------------------------------------------------------
! Returns fuel material effective creep rate (1/s)
!--------------------------------------------------------------------------------------------------
   real(c_double) function fcreep(i,j) bind(C,name='fcreep_')

   integer(c_int) i,j
   real(c_double) rate,sg

   fcreep=0.0d0
   if(fmat .eq. 'mox' .and. sigf(i,j) .gt. 0.0d0)then
      sg = sigf(i,j)*1.0d6
!     fission rate per unit volume (fiss/m3-s): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
      rate = qqv1(j)/(1.60214d-13*210.0d0)
      fcreep = fcreepc(1) * sg**fcreepc(2) * rate * dexp(-fcreepc(3)/(rmu*tem(i,j)))
   end if
   return
   end function fcreep

!--------------------------------------------------------------------------------------------------
! Calculation of the fuel Young's modulus, unit: MPa
!     input arguments :
!        tk .... temperature (k)
!--------------------------------------------------------------------------------------------------
   real(c_double) function felmod(tk,den,cont) bind(C,name='felmod_')

   real(c_double) tk,den,cont,tden,por

   if(fmat.eq.'mox')then
      tden=11460.d0*cont+10960.d0*(1.d0-cont)
      por=1.d0-den/tden
!     MATPRO correlation
      felmod=(2.334d11*(1.0d0-2.752d0*por)*(1.d0-1.0915d-4*tk))*(1.d0+0.15d0*cont)
   else
      write(*,*)'wrong fuel material:',fmat
      stop
   end if
   felmod = felmod/1.0d6
   return
   end function felmod

!--------------------------------------------------------------------------------------------------
! Calculation of fission gas release
!--------------------------------------------------------------------------------------------------
   real(c_double) function fgr(j)

   integer(c_int) j,i
   real(c_double) a,fgrR,fgrU,aR,r

!  Waltar and Reynolds model
   if(bup(j) .eq. 0.0d0)then
      fgr = 0.0d0
   else
!     calculate release fraction for the restructured fuel
      a = 4.7d0/bup(j) * (1.0d0-dexp(-bup(j)/5.9d0))
      fgrR = max(0.0d0, 1.0d0 - a)
   
!     calculate release fraction for the unrestructured fuel
      if(bup(j) .le. 3.5d0)then
          fgrU = 0.0d0
      else    
          a = 25.6d0/(bup(j)-3.5d0) * (1.0d0-dexp(-(bup(j)/3.5d0-1.0d0)))*dexp(-0.0125d0*(ql(j)/1.d3))
          if(bup(j) .ge. 49.2d0) a = a*exp(-0.3d0*(bup(j)-49.2d0))
          fgrU = max(0.0, 1.0d0 - a)
      end if
!     fractional fuel areas associated with restructured zone
      aR = 0.0d0
!     assume that the temperature 1500 C defines the boundary between unrestructured and restructured zones (Waltar & Reynolds, p.198).
      do i=1,nf
         if(tem(i,j)-273.15d0 .ge. 1500.0d0)aR = aR + az0(i)
      end do
      aR = aR/azf0
      fgr = (1.0d0 - aR)*fgrU + aR*fgrR
   end if
   
   return
   end function fgr

!--------------------------------------------------------------------------------------------------
! Calculation of the fuel local thermal conductivity, unit : w/(m*k)
!     input arguments :
!        dens - density (kg/m3)
!        bup - burn-up (mwd/kgu)
!        tk - temperature (k)
!        pucon - Pu fraction (-)
!        sto - fuel stoichiometry
!--------------------------------------------------------------------------------------------------
      real(c_double) function flamb(dens,bup,tk,pucon,sto) bind(C,name='flamb_')

      real(c_double) dens,bup,tk,pucon,sto,ac,tden,por

      if(fmat.eq.'mox')then
!        stoichometry
         ac=1.320d0*((2.0d0- sto)+0.0093d0)**0.5d0 - 0.091d0 + 0.0038d0*bup/9.33d0
         tden=11460.d0*pucon+10960.d0*(1.d0-pucon)
!        fuel porosity
         por=1.d0-dens/tden
!        FBR MOX fuel  thermal conductivity in [W/mK] according to Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
         flamb=(1.0d0/(ac + 2.493d-4*tk) + 88.4d-12*tk**3)*(1.d0-por)/(1.0d0+2.0d0*por)
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end function flamb

!--------------------------------------------------------------------------------------------------
! Calculation of the fuel Poisson ratio, unitless
!--------------------------------------------------------------------------------------------------
      real(c_double) function fpoir() bind(C,name='fpoir_')

!     MATPRO model
      if(fmat.eq.'mox')then
         fpoir=0.276d0
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end function fpoir

!--------------------------------------------------------------------------------------------------
! Volumetric fuel swelling rate (-)
!--------------------------------------------------------------------------------------------------
   real(c_double) function fswel(i,j) bind(C,name='fswel_')

   integer(c_int) i,j
   real(c_double) rate, dT, E, B

   fswel=0.0d0
   if(ifswel .eq. 1)then
      if(fmat .eq. 'mox')then
!        !fission rate per unit volume (fiss/m3-s): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
         rate = qqv1(j)/(1.60214d-13*210.0d0)
!        MATPRO model for swelling due to gaseous FPs
         if(tem(i,j) .lt. 2800.0d0)then
            dT = 2800.0d0-tem(i,j)
         !if(tem(i,j) .lt. 3700.0d0)then
         !   dT = 3700.0d0-tem(i,j)
!           total energy generated per unit volume (J/m3)
            E = bup(j)*1.0d6*8.64d4*rof0
!           total number of fissions per unit volume (fiss/m3): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
            B = E/(1.60214d-13*210.0d0)
            fswel = 8.8d-56 * dT**11.73d0 * exp(-0.0162d0*dT) * exp(-8.0d-27*B) * rate
!           correction to account for reduction of swelling with increase of fission gas release
            !fswel = fswel / (1.d0 + dexp(50.d0*(fgrj(j) - 0.9d0)))
         else
            fswel = 0.0d0
         end if

!        MATPRO model for swelling due to solid FPs
         fswel = fswel + 2.5d-29*rate
      end if
      if(fswelmlt .ne. 0.0d0)fswel = fswel*fswelmlt
   end if
   return
   end function fswel

!--------------------------------------------------------------------------------------------------
! Calculation of fuel thermal expansion, unit : m/m
!       input arguments :
!          tk .... temperature (k)
!--------------------------------------------------------------------------------------------------
   real(c_double) function ftexp(tk,cont) bind(C,name='ftexp_')

   real(c_double) tk,cont,tc,ftexpUO2,ftexpPuO2

   if(fmat.eq.'mox')then
!     MATPRO model
      ftexpUO2  = 1.0d-5*tk - 3.0d-3 + 4.d-2*dexp(-5000.d0/tk)
      ftexpPuO2 = 0.9d-5*tk - 2.7d-3 + 7.d-2*dexp(-5072.d0/tk)
      ftexp = ftexpUO2*(1.0-cont) + ftexpPuO2*cont 
   else
      write(*,*)'wrong fuel material:',fmat
      stop
   end if
   return
   end function ftexp

!--------------------------------------------------------------------------------------------------
! Open and closed fuel/clad gap conductance
!  j - axial slice
!--------------------------------------------------------------------------------------------------
   real(c_double) function gaphtc(j) bind(C,name='gaphtc_')

   integer(c_int) j,i,jj
   real(c_double) flamb,clamb
   real(c_double) gtemp,fe,sbc,emissf,emissc,cmhard,r,conf,conc,fkm,ahe,pgas,gask_(4), &
          xmol_(4),sumx,tf,tc,hgap(3),mw(4),gap_
!          he        ar        kr      xe
   data mw/4.0026d0, 39.948d0, 83.8d0, 131.30d0/

!  Stefan-Boltzmann constant (w/m**2-k**4)
   data sbc/5.6697d-8/

   tf=tem(nf,j)
   tc=tem(nf+1,j)

!  average gap temperature
   gtemp=(tf+tc)/2.d0
!  helium thermal conductivity
   gask_(1)=2.639d-3*gtemp**0.7085d0
!  argon thermal conductivity
   gask_(2)=2.986d-4*gtemp**0.7224d0
!  kripton thermal conductivity
   gask_(3)=8.247d-5*gtemp**0.8363d0
!  xenon thermal conductivity
   gask_(4)=4.351d-5*gtemp**0.8618d0
   if(gmat.eq.'he')then
      xmol_(1)=mu0
      xmol_(2)=0.d0
   else ! if(gmat.eq.'ar')then
      xmol_(1)=0.d0
      xmol_(2)=mu0
   end if

!  number of moles of released fission gases (0.8846xe,0.0769kr,0.0385he)
   xmol_(1) = xmol_(1) + 0.0385d0*fgrel
   xmol_(3) = 0.0769d0*fgrel
   xmol_(4) = 0.8846d0*fgrel

   sumx=0.d0
   do i=1,4
      sumx=sumx+xmol_(i)
      xmol(i)=xmol_(i)
   end do
   gask(j)=1.0d0
   do i=1,4
      gask(j) = gask(j)*gask_(i)**(xmol_(i)/sumx)
   end do

!  open gap conductance
   ahe=0.425d0-2.3d-4*dmin1(gtemp,1000.d0)
   pgas=dmax1(gpres,0.1d0)*1.0d6
!  accomodation distance
   ajump(j)=0.024688d0*gask(j)*dsqrt(gtemp)*2.d0/(pgas*ahe)

!  using FRAPCON model reduce the gap width due to fuel cracking while calculating the gap conductance 
   if(ifreloc .eq. 1)then
      reloc(j)=0.25d0 + dmax1( 5.0d0, dmin1(10.d0,ql(j)/4000.d0) )*( dmin1( 1.0d0, bup(j)/5.d0) + 1.0d0 )/100.d0
   else
      reloc(j)=0.0d0
   end if

!  effective roughness
   r=sqrt(rufc**2+ruff**2)
   gap_=max(gap(j), r)
 
   gapth(j) = gap_*(1.0d0 - reloc(j)) + ajump(j)
   hgap(1)=gask(j)/gapth(j)

!  radiation conductance
!  set fuel and clad surface emissivities to reasonable value.
!  according to KIT recommendation (Struwe+Schikorr) 15.06.2011 
   emissf = 0.8d0
   emissc = 0.9d0

   fe=1.0d0/(1.0d0/emissf+(rfo0/rci0)*(1.0d0/emissc-1.0d0))
   hgap(2)=sbc*fe*(tf**2+tc**2)*(tf+tc)

   if(flag(j).eq.'open')then
      hgap(3)=0.d0
   else
!     average cladding temperature
      tc = 0.0d0
      do i=1,nc
         tc = tc + tem(nf+i,j)*az0(nf+i)
      end do
      tc = tc/azc0
!     Meyer hardness for SS (from SAS4A hgap routine)
      if(tc.le.893.9203d0)then
         cmhard=5.961d9*tc**(-0.206d0)
      else
         cmhard=2.750d28*tc**(-6.530d0)
      end if
      if(cmhard.lt.1.d5)cmhard=1.d5

!     effective roughness
      r=sqrt(rufc**2+ruff**2)

!     effective conductivity
      conf=flamb(rof0,bup(j),tf,pucont,sto0)
      conc=clamb(tc)
      fkm=2.d0*conf*conc/(conf+conc)

!     closed gap conductance
      hgap(3)=33.3d0*fkm*(pfc(j)-gpres)*1.0d6/(cmhard*sqrt(r))
   end if
   gaphtc = hgap(1) + hgap(2) + hgap(3)

   hgapi(1,j)=hgap(1)
   hgapi(2,j)=hgap(2)
   hgapi(3,j)=hgap(3)
   
   return
   end function gaphtc

!--------------------------------------------------------------------------------------------------
! Internal gas pressure calculation
!--------------------------------------------------------------------------------------------------
   real(c_double) function gpresf(time) bind(C,name='gpresf_')

   integer(c_int) j
   real(c_double) polp
   real(c_double) time,v_t,tgap,thol,tp

!  table of gas plenum temperature (K) vs time
   if(ntple.gt.0)then
      tp = polp(time,tple,ttple,ntple)
   else
      tp = tem0
   end if

!  volume-to-temperature ratio
   v_t = vgp/tp
   do j=1,nzz
!    gap and hole temperature
     tgap = (tem(nf,j) + tem(nf+1,j))/2.d0
     thol = tem(1,j)
!    gap vol/tem + hole vol/tem
     v_t = v_t + pi*dz0(j)*((rfo0**2 - rci0**2)/tgap + rfi0**2/thol)
   end do

!  free volume gas pressure (MPa)
   gpresf = (mu0 + fgrel)*rmu / v_t * 1d-6

   return
   end function gpresf

!-----------------------------------------------------------------------------------------------
! Write outfrd file
!-----------------------------------------------------------------------------------------------
   subroutine outfrd(time,nnn)

   integer(c_int) nnn
   integer(c_int) i,j
   real(c_double) time
   real(c_double) flamb
   character*12 s
   character*20 title

   write(s,'(i12.12)') nnn
   open(700,file='outfrd' // s)
   title = "time (s)            "
   write(700,12)title,time
   
   title = "time (d)            "
   write(700,12)title,time/86400.d0
   
   write(700,'(500a)')('-',j=1,100)
   title = "rof0 (kg/m3)        "
   write(700,11)title,rof0
   
   title = "fggen (cm3)         "
   write(700,11)title,fggen*rmu*293.15d0/1.0d5 * 1.0d6
   
   title = "fgrel (cm3)         "
   write(700,11)title,fgrel*rmu*293.15d0/1.0d5 * 1.0d6
   
   title = "fgrel (%)           "
   write(700,11)title,fgrel/(fggen+1.d-21)*100.d0
   
   title = "gpres (MPa)         "
   write(700,11)title,gpres
   
   write(700,'(500a)')('-',j=1,100)
   title = "z (m)               "
   write(700,11)title,(dz0(j)*(real(j)-0.5d0),j=1,nzz)
   
   title = "tfin (C)            "
   write(700,11)title,(tem(1,j)-273.15d0,j=1,nzz)
   
   title = "tfout (C)           "
   write(700,11)title,(tem(nf,j)-273.15d0,j=1,nzz)
   
   title = "tcin (C)            "
   write(700,11)title,(tem(nf+1,j)-273.15d0,j=1,nzz)
   
   title = "tcout (C)           "
   write(700,11)title,(tem(nf+nc,j)-273.15d0,j=1,nzz)
   
   title = "rfi (m)             "
   write(700,11)title,(rfi(j),j=1,nzz)
   
   title = "rfo (m)             "
   write(700,11)title,(rfo(j),j=1,nzz)
   
   title = "rci (m)             "
   write(700,11)title,(rci(j),j=1,nzz)
   
   title = "rco (m)             "
   write(700,11)title,(rco(j),j=1,nzz)
   
   title = "dzf (m)             "
   write(700,11)title,(dzf(j),j=1,nzz)
   
   title = "dzc (m)             "
   write(700,11)title,(dzc(j),j=1,nzz)
   
   title = "reloc (-)           "
   write(700,11)title,(reloc(j),j=1,nzz)
   
   title = "qv (W/m3)           "
   write(700,11)title,(qqv1(j),j=1,nzz)
   
   title = "ql (W/cm)           "
   write(700,11)title,(ql(j)/100.d0,j=1,nzz)
   
   title = "ql2 (W/cm)          "
   write(700,11)title,(ql(j)/100.d0*dz0(j)/dzf(j),j=1,nzz)

   title = "bup (MWd/kg)        "
   write(700,11)title,(bup(j),j=1,nzz)

   title = "gap (m)             "
   write(700,11)title,(gap(j),j=1,nzz)
   
   title = "gapth (m)             "
   write(700,11)title,(gapth(j),j=1,nzz)
   
   title = "hgap (W/m2K)        "
   write(700,11)title,(hgapt(j),j=1,nzz)
   
   title = "hgap1 (W/m2K)       "
   write(700,11)title,(hgapi(1,j),j=1,nzz)
   
   title = "hgap2 (W/m2K)       "
   write(700,11)title,(hgapi(2,j),j=1,nzz)
   
   title = "hgap3 (W/m2K)       "
   write(700,11)title,(hgapi(3,j),j=1,nzz)
   
   title = "ajump (m)           "
   write(700,11)title,(ajump(j),j=1,nzz)
   
   title = "gask (W/m-K)        "
   write(700,11)title,(gask(j),j=1,nzz)
   
   title = "pfc (MPa)           "
   write(700,11)title,(pfc(j),j=1,nzz)
   
   write(700,'(500a)')('-',j=1,100)
   title = "r (m)               "
   write(700,11)title,(rad0(i),i=1,nf+nc)
   
   title = "temperature (C)     "
   do j=1,nzz
      write(700,10)title,j,(tem(i,j)-273.15d0,i=1,nf+nc)
   end do
   
   title = "sig h (MPa)         "
   do j=1,nzz
      write(700,10)title,j,(sigfh(i,j),i=1,nf),(sigh(i,j),i=1,nc)
   end do
   
   title = "sig r (MPa)         "
   do j=1,nzz
      write(700,10)title,j,(sigfr(i,j),i=1,nf),(sigr(i,j),i=1,nc)
   end do
   
   title = "sig z (MPa)         "
   do j=1,nzz
      write(700,10)title,j,(sigfz(i,j),i=1,nf),(sigz(i,j),i=1,nc)
   end do
   
   title = "eps h total (%)     "
   do j=1,nzz
      write(700,10)title,j,(efh(i,j)*100.d0,i=1,nf),(eh(i,j)*100.d0,i=1,nc)
   end do
   
   title = "eps r total (%)     "
   do j=1,nzz
      write(700,10)title,j,(efr(i,j)*100.d0,i=1,nf),(er(i,j)*100.d0,i=1,nc)
   end do
   
   title = "eps z total (%)    "
   do j=1,nzz
      write(700,10)title,j,efz(j)*100.d0,ez(j)*100.d0
   end do

   title = "eps h creep (%)     "
   do j=1,nzz
      write(700,10)title,j,(efch(i,j)*100.0d0,i=1,nf),(ech(i,j)*100.0d0,i=1,nc)
   end do

   title = "eps r creep (%)     "
   do j=1,nzz
      write(700,10)title,j,(efcr(i,j)*100.0d0,i=1,nf),(ecr(i,j)*100.0d0,i=1,nc)
   end do
   
   title = "eps z creep (%)     "
   do j=1,nzz
      write(700,10)title,j,(efcz(i,j)*100.0d0,i=1,nf),(ecz(i,j)*100.0d0,i=1,nc)
   end do
   
   title = "eps thermal lin (%) "
   do j=1,nzz
      write(700,10)title,j,(eft(i,j)*100.0d0,i=1,nf),(et(i,j)*100.0d0,i=1,nc)
   end do
   
   title = "eps swell vol (%)   "
   do j=1,nzz
      write(700,10)title,j,(efs(i,j)*100.0d0,i=1,nf)
   end do

   close(700)
10 format(a20,1x,"iz:",i2,500(1x,1pe12.5))
11 format(a20,1x,"   ",2x,500(1x,1pe12.5))
12 format(a20,1x,"   ",2x,500(1x,1pe12.5))

   if(time .lt. 60.d0)then
      write(*,'("Time: ",(1pe12.5)," s  ")')time
   else if(time .lt. 60.d0*60.d0)then
      write(*,'("Time: ",(1pe12.5)," min")')time/60.d0
   else if(time .lt. 60.d0*60.d0*24.d0)then
      write(*,'("Time: ",(1pe12.5)," h  ")')time/3600.d0
   else
      write(*,'("Time: ",(1pe12.5)," d  ")')time/(3600.d0*24.d0)
   end if
   return
   end


!--------------------------------------------------------------------------------------------------
! Write or read rstfrd file
!--------------------------------------------------------------------------------------------------
   subroutine rstfrd(act,time,nnn,y,yp)

   integer(c_int), intent(inout) :: nnn
   real(c_double), intent(inout) :: time, y(maxeq), yp(maxeq)
   integer(c_int) i, ios
   character*5 act
   character*12 s
   character*20 title
   logical isfile


   if(act .eq. 'write')then
      write(s,'(i12.12)') nnn
      open(700,file='rstfrd' // s)

      title = "time (s)            "
      write(700,'(a8,1x,10000(1pe16.9,1x))')title,time
      
      title = "time (d)            "
      write(700,'(a8,1x,10000(1pe16.9,1x))')title,time/86400.d0
      
      title = "y                   "
      write(700,'(a8,1x,10000(1pe16.9,1x))')title,(y(i),i=1,neq)
      
      title = "yp                  "
      write(700,'(a8,1x,10000(1pe16.9,1x))')title,(yp(i),i=1,neq)
   else if(act .eq. 'read ')then
      write(s,'(i12.12)') nnn
!     check if the rstfrd file exists
      inquire(file='rstfrd' // s, exist=isfile)
      if(.not. isfile)then
         write(*,*)'ERROR: restart file does not exist ', 'rstfrd' // s
         stop
      end if
      open(700,file='rstfrd' // s)
      do
         read(700,'(a8)', iostat=ios)title
         if (ios /= 0) exit
         if(title(:8) .eq. 'time (s)')then
            backspace 700
            read(700,'(a8,1x,1pe16.9)')title,time
         else if(title(:2) .eq. 'y ')then
            backspace 700
            read(700,'(a8,1x,10000(1pe16.9,1x))')title,(y(i),i=1,neq)
         else if(title(:2) .eq. 'yp')then
            backspace 700
            read(700,'(a8,1x,10000(1pe16.9,1x))')title,(yp(i),i=1,neq)
         end if
      end do
   end if

   close(700)
   return
   end

!--------------------------------------------------------------------------------------------------
! polp is a linear interpolation function whose value is equal to
! interpolated value.
! tc-is the value of the independent variable
! crf-dependent variable array
! crft-independent variable array
! n-number of points in variable arrays
!-------------------------------------------------------------------------------------------------
   real(c_double) function polp(tc,crf,crft,n) bind(C,name='polp_')

   integer(c_int) n,i,j
   real(c_double) tc,crf(n),crft(n)
   i=1
   if(tc.ge.crft(1))then
      i=n
      if(tc.lt.crft(n))then
         do j=1,n
            i=j
            if(tc-crft(j) < 0.0d0)then
               polp=crf(i-1)+(crf(i)-crf(i-1))*(tc-crft(i-1))/(crft(i)-crft(i-1))
               return
            else if(tc-crft(j) == 0.0d0)then
               polp=crf(i)
               return
            end if
         end do
      end if
   end if
   polp=crf(i)
   return
   end function polp

   real(c_double) function polp2(xa,x,y,n,eps) result(ya) bind(C,name='polp2_')
   integer(c_int) n,i,j
   real(c_double) xa,x(maxtab),y(maxtab),eps,dx,dydx,dydx_,alfa

   if(xa .lt. x(1))then
      ya = y(1)
   else if(xa .gt. x(n))then
      ya = y(n)
   else
!     there are n (x,y) couples and n-1 intervals
      i = 1
!     find the interval i in which xa is
      do
         if(xa .lt. x(i)) exit
         i = i + 1
      end do
      dydx = (y(i) - y(i-1))/(x(i) - x(i-1))
      ya = y(i-1) + dydx*(xa-x(i-1))
      if(xa .lt. x(i-1)+eps)then
         alfa = (xa - (x(i-1) - eps))/(2.0d0*eps)
         if(i-1 .gt. 1)then
            dydx_ = (y(i-1) - y(i-2))/(x(i-1) - x(i-2))
            ya = alfa*ya + (1.0d0 - alfa)*(y(i-2) + dydx_*(xa-x(i-1)))
         end if
      else if(xa .gt. x(i)-eps)then
         alfa = (xa - (x(i) - eps))/(2.0d0*eps)
         if(i .eq. n)then
            dydx_ = 0.0d0
         else
            dydx_ = (y(i+1) - y(i))/(x(i+1) - x(i))
         end if
         ya = (1.0d0 - alfa)*ya + alfa*(y(i+1) + dydx_*(xa-x(i-1)))
      end if
   end if

   end function polp2

!--------------------------------------------------------------------------------------------------
! Reads input deck and assigns initial values
!--------------------------------------------------------------------------------------------------
      subroutine rdinp()

      integer(c_int) numlin,i1,i2,i,j
      real(c_double) ftexp,felmod_,fpoir,ctexp,celmod,cpoir,gaphtc_
      real(c_double) r1,eft0,fyng0,fpoir0,ef0,sigf0,et0,cyng0,cpoir0,tfree0,vfree0
      character*1 c
      character*20 w2

!     table dimensions
      nqv=0
      ntco=0
      nfsto=0
      nzz=0
!     default values for options      
      ifgr=0
      ifcreep=0
      ifswel=0
      iccreep=0
      ifreloc=0
      inomech=0
      irst=-1

      open(700, file='fred.inp')
!     first symbol of the card
      read(700,'(a1)',err=1000) c
!     input deck line number
      numlin = 1

!     $ - terminal card
      do while(c .ne. '$')
!        * - comment card
         if(c .ne. '*') then
            backspace 700
!           first integer of the card
            read(700,*,err=1000)i1

!           options card
            if(i1.eq.000001) then
               backspace 700
               read(700,*,err=1000)i1,w2
               
               if(w2 .eq. 'FGR')then
                  ifgr=1
                  
               else if(w2 .eq. 'FUEL_CREEP')then
                  ifcreep=1
                  
               else if(w2 .eq. 'FUEL_SWEL')then
                  ifswel=1
                  backspace 700
                  read(700,*,err=1000)i1,w2,fswelmlt

               else if(w2 .eq. 'CLAD_CREEP')then
                  iccreep=1

               else if(w2 .eq. 'FUEL_RELOC')then
                  ifreloc=1

               else if(w2 .eq. 'NOMECH')then
                  inomech=1

               else if(w2 .eq. 'RESTART')then
                  backspace 700
                  read(700,*,err=1000)i1,w2,irst

               end if   

!           initial temperature card
            else if(i1.eq.000002) then
               backspace 700
               read(700,*,err=1000)i1,tem0

!           fixed time card
            else if(i1.eq.000004) then
               backspace 700
               ntfix = ntfix + 1
               read(700,*,err=1000)i1,tfix(ntfix)

!           fuel rod axial division card
            else if(i1.eq.000003)then
               backspace 700
!              dz0 - axial slice height              
               read(700,*,err=1000)i1,r1,i2
               if(i2.gt.maxz)stop 'card 000003:increase maxz'
               nzz=i2
               do j=1,i2
                  dz0(j)=r1
               end do

!           fuel card
            else if(i1.eq.100001)then
               backspace 700
!              fmat  - fuel material (MOX)
!              rof - fuel density (kg/m3)
!              pucont - plutonium content
!              rfi0  - inner fuel radius (m)
!              rfo0  - outer fuel radius (m)
!              ruff - arithmetic mean roughness height of fuel (m)
!              sto0 - initial fuel stoichiometry
!              nf    - number of radial fuel nodes
               read(700,*,err=1000) i1,fmat,rof0,pucont,rfi0,rfo0,ruff,sto0,nf
               if(fmat.eq.'none')nf=0
               do j=1,nzz
                  do i=1,nf
                     rof(i,j)=rof0
                  end do
               end do

 !          gap card
            else if(i1.eq.100002)then
               backspace 700
!              gmat - gap material (He, Ar)
!              gap  - gap width (m)
!              gpres0 - initial inner gas pressure (MPa)
!              vgp - gas plenum volume
               read(700,*,err=1000)i1,gmat,gap(1),gpres0,vgp
               if(gmat.ne.'he' .and. gmat.ne.'ar')then
                  write(*,*)'wrong gas gap material ',gmat
                  write(*,*)'only he and ar allowed'
                  stop
               end if
               do j=2,nzz
                  gap(j)=gap(1)
               end do

!           clad card
            else if(i1.eq.100003)then
               backspace 700
!              cmat - clad material (aim1, t91)
!              rco0  - outer clad radius (m)
!              roc - clad density (kg/m3)
!              rufc - arithmetic mean roughness height of cladding (m)
!              nc - number of clad nodes
               read(700,*,err=1000)i1,cmat,rco0,roc0,rufc,nc

!           time card
            else if(i1.eq.200000)then
               backspace 700
!              tend - end time
!              dtout - output time step
!              hmax - maximum time step
               read(700,*,err=1000)i1,tend,dtout,hmax

!           clad outer temperature card
            else if(i1.eq.200001)then
               ntco=ntco+1
               if(ntco.gt.maxtab) stop 'increase maxtab'
               backspace 700
!              ntco  - number of clad outer temperature table entries
!              tco - array of fuel rod clad outer temperature vs time ttco (K)
!              ttco - time array corresponding to tco (s)
               read(700,*,err=1000)i1,ttco(ntco),(tco(i,ntco),i=1,nzz)
               
!           power density card
            else if(i1.eq.200002)then
               nqv=nqv+1
               if(nqv.gt.maxtab) stop 'increase maxtab'
               backspace 700
!              nqv  - number of power table entries
!              qv - array of fuel rod power density vs time tqv (kW/m)
!              tqv - time array corresponding to qv (s)
               read(700,*,err=1000)i1,tqv(nqv),(qv(i,nqv),i=1,nzz)

!           fuel creep rate law constants
            else if(i1.eq.200003) then
               backspace 700
!              ECR = fcreepc(1) * SIG**fcreepc(2) * exp(-fcreepc(3)/(R*T))
!              where ECR is the fuel creep rate (1/s), 
!                    SIG the effective (Von Mises) stress (Pa), 
!                    R the universal gas constant (J/mol-K),
!                    T the temperature (K)
               read(700,*,err=1000)i1,fcreepc(1),fcreepc(2),fcreepc(3)

!           fuel stoichiometry vs fuel burnup table
            else if(i1.eq.200004) then
               nfsto=nfsto+1
               if(nfsto.gt.maxtab) stop 'increase maxtab'
               backspace 700
!              nfsto - number of fuel stoichiometry vs fuel burnup entries
!              bsto - burnup array corresponding to stob
!              stob - array of fuel stoichiometry vs fuel burnup
               read(700,*,err=1000)i1,bsto(nfsto),stob(nfsto)

!           coolant pressure card
            else if(i1.eq.200005)then
               backspace 700
!              pcool - coolant pressure (MPa)
               read(700,*,err=1000)i1,pcool

!           gas plenum temperature card
            else if(i1.eq.200006)then
               ntple=ntple+1
               if(ntple.gt.maxtab) stop 'increase maxtab'
               backspace 700
!              tple - gas plenum temperature (K)
               read(700,*,err=1000)i1,ttple(ntple),tple(ntple)

            else
               write(*,*)'input error: wrong card number',i1
               stop
            end if
         end if
!        first symbol of the next string
         read(700,'(a1)',err=1000)c
         numlin=numlin+1
      end do
      close(700)

!     some preparations...

!     outer cladding temperature table
      if(ntco .eq. 0)then
         stop 'INPUT ERROR: Specify cladding temperature (card 200000)'
      end if

!     set initial geometry
      do j=1,nzz
         rci0=rfo0+gap(j)
         rfi(j)=rfi0
         rfo(j)=rfo0
         rco(j)=rco0
         rci(j)=rci0
         dzf(j)=dz0(j)
         dzc(j)=dz0(j)
         flag(j)='open'
      end do

      drf0=(rfo0-rfi0)/dfloat(nf-1)
!     radii of fuel nodes
      do i=1,nf
         if(i.eq.1)then
            rad0(1)=rfi0
         else
            rad0(i)=rad0(i-1)+drf0
         end if
         do j=1,nzz
            drf(i,j)=drf0
         end do
      end do

!     radii of clad nodes
      drc0=(rco0-rci0)/dfloat(nc-1)
      do i=nf+1,nf+nc
         if(i.eq.nf+1)then
            rad0(nf+1)=rci0
         else
            rad0(i)=rad0(i-1)+drc0
         end if
      end do
      do j=1,nzz
      do i=1,nf+nc
         rad(i,j)=rad0(i)
      end do
      end do
      do j=1,nzz
      do i=1,nc-1
         drc(i,j)=drc0
      end do
      end do

!     initial radii of boundaries between nodes
      do i=1,nf+nc-1
         rad_0(i)=(rad0(i) + rad0(i+1))/2.0d0
      end do

      do i=1,nf
!        XS area
         if(i.eq.1) then
            az0(1)=pi*((rfi0+drf0/2)**2-rfi0**2)
         else if(i.eq.nf) then
            az0(nf)=pi*(rfo0**2-(rfo0-drf0/2.0d0)**2)
         else
            az0(i)=pi*(rad0(i)+drf0/2.0d0)**2 - pi*(rad0(i)-drf0/2.0d0)**2
         end if
      end do

!     initial fuel stack height
      h0 = 0.d0
      do j=1,nzz
         h0 = h0 + dz0(j)
      end do
!     fuel pellet cross section area
      azf0=pi*(rfo0**2-rfi0**2)

      do i=nf+1,nf+nc
!        axial area of volume
         if(i.eq.nf+1) then
            az0(i)=pi*((rci0+drc0/2.0d0)**2-rci0**2)
         else if(i.eq.nf+nc) then
            az0(i)=pi*(rco0**2-(rco0-drc0/2.0d0)**2)
         else
            az0(i)=pi*(rad0(i)+drc0/2.0d0)**2 - pi*(rad0(i)-drc0/2.0d0)**2
         end if
      end do
!     clad cross section area
      azc0=pi*(rco0**2-rci0**2)

      fggen=0.d0
      fgrel=0.d0

!     set initial temperatures
      do j=1,nzz
         do i=1,nf+nc
            tem(i,j)=tem0
         end do
      end do

!     as-fabricated free volume
      vfree0 = pi*(rci0**2 - rfo0**2 + rfi0**2)*h0 + vgp
!     as-fabricated temperature
      tfree0 = 293.d0
!     as fabricated gas amount
      mu0 = gpres0*1.0d6*vfree0/tfree0/rmu

!     initial gas pressure
      gpres = gpresf(0.0d0)
!     set initial fuel stresses
      sigf0 = -gpres
      do j=1,nzz
         do i=1,nf
            sigfh(i,j) = sigf0
            sigfr(i,j) = sigf0
            sigfz(i,j) = sigf0
         end do
      end do
!     set initial fuel strains
      eft0 = ftexp(tem0,pucont) - ftexp(293.15d0,pucont)
      fyng0 = felmod(tem0,rof0,pucont)
      fpoir0 = fpoir()
      ef0 = eft0 + sigf0*(1.0d0 - 2.0d0*fpoir0) / fyng0
      do j=1,nzz
         efz(j)=ef0
         do i=1,nf
            efh(i,j) = ef0
            efr(i,j) = ef0
            eft(i,j) = eft0
         end do
      end do
!     set initial clad stresses
      do j=1,nzz
         do i=1,nc
            sigh(i,j) = (gpres - pcool)*(1.0d0 + rci0**2/rad(nf+i,j)**2)/(1.0d0 - rci0**2/rco0**2) - gpres
            sigr(i,j) = (gpres - pcool)*(1.0d0 - rci0**2/rad(nf+i,j)**2)/(1.0d0 - rci0**2/rco0**2) - gpres
            sigz(i,j) = (gpres*rci0**2 - pcool*rco0**2)/(rco0**2 - rci0**2)
         end do
      end do
!     set initial clad strains
      et0 = ctexp(tem0) - ctexp(293.15d0)
      cyng0 = celmod(tem0)
      cpoir0 = cpoir()
      do j=1,nzz
         ez(j) = et0 + (sigz(1,j) - cpoir0*(sigr(1,j)+sigh(1,j)))/cyng0
         do i=1,nc
            eh(i,j) = et0 + (sigh(i,j) - cpoir0*(sigr(i,j)+sigz(i,j)))/cyng0
            er(i,j) = et0 + (sigr(i,j) - cpoir0*(sigh(i,j)+sigz(i,j)))/cyng0
            et(i,j) = et0
         end do
      end do
!     update geometry and set initial gap conductance
      do j=1,nzz
         call update_geom(j)
         gap(j) = rci(j) - rfo(j)
         hgapt(j) = gaphtc(j)
      end do

!     number of equations
      neq = 0
      neq = neq + 1                  !fggen
      neq = neq + 1                  !fgrel
      neq = neq + 1                  !gpres

      neq_j = 0
      neq_j = neq_j + 1              !power density
      neq_j = neq_j + 1              !burnup
      neq_j = neq_j + 1              !gap
      neq_j = neq_j + 1              !pfc
      neq_j = neq_j + 1              !hgapt
      neq_j = neq_j + nf             !fuel temperature
      neq_j = neq_j + nc             !clad temperature
      if(inomech .eq. 0) then
         neq_j = neq_j + nf          !efs
         neq_j = neq_j + nf          !eft
         if(ifcreep .eq. 1)then
            neq_j = neq_j + nf       !efce
            neq_j = neq_j + nf       !efch
            neq_j = neq_j + nf       !efcr
            neq_j = neq_j + nf       !efcz
         end if
         neq_j = neq_j + nf          !efh
         neq_j = neq_j + nf          !efr
         neq_j = neq_j + 1           !efz
         neq_j = neq_j + nf          !sigfh
         neq_j = neq_j + nf          !sigfr
         neq_j = neq_j + nf          !sigfz
         neq_j = neq_j + nc          !et
         if(iccreep .eq. 1)then
            neq_j = neq_j + nc       !ece
            neq_j = neq_j + nc       !ech
            neq_j = neq_j + nc       !ecr
            neq_j = neq_j + nc       !ecz
         end if
         neq_j = neq_j + nc          !eh
         neq_j = neq_j + nc          !er
         neq_j = neq_j + 1           !ez
         neq_j = neq_j + nc          !sigh
         neq_j = neq_j + nc          !sigr
         neq_j = neq_j + nc          !sigz
      end if

      neq = neq + neq_j*nzz

      if(neq .gt. maxeq)then
         write(*,*)'increase maxeq in global.f90 up to ', neq
         stop
      end if

      return
!     input error 
1000  write(*,*)'input error: line number: ', numlin
      stop
      end

!--------------------------------------------------------------------------------------------------
! Update geometry for axial layer j
!--------------------------------------------------------------------------------------------------
   subroutine update_geom(j)

   integer(c_int) i,j,k
   real(c_double) efr_(maxr),er_(maxr)

   do i=1,nf-1
      efr_(i) = 0.5d0*(efr(i,j) + efr(i+1,j))
   end do
   do i=1,nc-1
      er_(i) = 0.5d0*(er(i,j) + er(i+1,j))
   end do

!  Given fuel strains (efh, efr, efz), calculate rad, drf, dzf
   do i=1,nf
!     fuel nodes radii
      rad(i,j) = rad0(i)*(1.d0 + efh(i,j))
      if(rad(i,j) .lt. 0.0d0)then
         write(*,*)rad(i,j), ': rad(i,j) < 0 (',i,j,')'
         stop
      end if
   end do
   do i=1,nf-1
!     fuel nodes thickness
      drf(i,j) = drf0*(1.0d0 + efr_(i))
      if(drf(i,j) .lt. 0.0d0)then
         write(*,*)drf(i,j), ': drf(i,j) < 0 (',i,j,')'
         stop
      end if
   end do
   dzf(j) = dz0(j)*(1.0d0 + efz(j))
   if(dzf(j) .lt. 0.0d0)then
      write(*,*)dzf(j), ': dzf(j) < 0 (',j,')'
      stop
   end if
   rfi(j) = rad(1,j)
   rfo(j) = rad(nf,j)
 
!  Given clad strains (eh, er, ez), calculate rad, drc, dzc
   do i=1,nc
!     clad nodes radii
      rad(nf+i,j) = rad0(nf+i)*(1.0d0 + eh(i,j))
      if(rad(nf+i,j).lt.0.0d0)then
         write(*,*)rad(nf+i,j), ': rad(nf+i,j) < 0 (',i,j,')'
         stop
      end if
   end do
   do i=1,nc-1
!     clad nodes thickness
      drc(i,j) = drc0*(1.0d0 + er_(i))
      if(drc(i,j).lt.0.0d0)then
         write(*,*)drc(i,j), ': drc(i,j) < 0 (',i,j,')'
         stop
      end if
   end do
   dzc(j) = dz0(j)*(1.0d0 + ez(j))
   if(dzc(j).lt.0.0d0)then
      write(*,*)dzc(j), ': dzc(j) < 0 (',j,')'
      stop
   end if
   rci(j) = rad(nf+1,j)
   rco(j) = rad(nf+nc,j)
   return
   end

!--------------------------------------------------------------------------------------------------
! Write variables to vector y
!--------------------------------------------------------------------------------------------------
   subroutine write_to_y(y)

   integer(c_int) i,j,k
   real(c_double), intent(inout) :: y(maxeq)

   k = 0
   k = k + 1
   y(k) = fggen
   k = k + 1
   y(k) = fgrel
   k = k + 1
   y(k) = gpres
   do j=1,nzz
      k = k + 1
      y(k) = qqv1(j)
      k = k + 1
      y(k) = bup(j)*8.64d10
      k = k + 1
      y(k) = gap(j)
      k = k + 1
      y(k) = pfc(j)
      k = k + 1
      y(k) = hgapt(j)
      do i=1,nf+nc
         k = k + 1
         y(k) = tem(i,j)
      end do
      if(inomech .eq. 0)then
         do i=1,nf
            k = k + 1
            y(k) = efs(i,j)
         end do
         do i=1,nf
            k = k + 1
            y(k) = eft(i,j)
         end do
         if(ifcreep .eq. 1)then
            do i=1,nf
               k = k + 1
               y(k) = efce(i,j)
            end do
            do i=1,nf
               k = k + 1
               y(k) = efch(i,j)
            end do
            do i=1,nf
               k = k + 1
               y(k) = efcr(i,j)
            end do
            do i=1,nf
               k = k + 1
               y(k) = efcz(i,j)
            end do
         end if
         do i=1,nf
            k = k + 1
            y(k) = efh(i,j)
         end do
         do i=1,nf
            k = k + 1
            y(k) = efr(i,j)
         end do
         k = k + 1
         y(k) = efz(j)
         do i=1,nf
            k = k + 1
            y(k) = sigfh(i,j)
         end do
         do i=1,nf
            k = k + 1
            y(k) = sigfr(i,j)
         end do
         do i=1,nf
            k = k + 1
            y(k) = sigfz(i,j)
         end do
         do i=1,nc
            k = k + 1
            y(k) = et(i,j)
         end do
         if(iccreep .eq. 1)then
            do i=1,nc
               k = k + 1
               y(k) = ece(i,j)
            end do
            do i=1,nc
               k = k + 1
               y(k) = ech(i,j)
            end do
            do i=1,nc
               k = k + 1
               y(k) = ecr(i,j)
            end do
            do i=1,nc
               k = k + 1
               y(k) = ecz(i,j)
            end do
         end if
         do i=1,nc
            k = k + 1
            y(k) = eh(i,j)
         end do
         do i=1,nc
            k = k + 1
            y(k) = er(i,j)
         end do
         k = k + 1
         y(k) = ez(j)
         do i=1,nc
            k = k + 1
            y(k) = sigh(i,j)
         end do
         do i=1,nc
            k = k + 1
            y(k) = sigr(i,j)
         end do
         do i=1,nc
            k = k + 1
            y(k) = sigz(i,j)
         end do
      end if
   end do
   return
   end

!--------------------------------------------------------------------------------------------------
! Read variables from vector y
!--------------------------------------------------------------------------------------------------
   subroutine read_from_y(y)

   integer(c_int) i,j,k
   real(c_double) y(maxeq)

   k = 0
   k = k + 1
   fggen = y(k)
   k = k + 1
   fgrel = max(0.0d0,y(k))
   k = k + 1
   gpres = y(k)
   do j=1,nzz
      k = k + 1
      qqv1(j) = y(k)
      k = k + 1
      bup(j) = y(k)/8.64d10
      k = k + 1
      gap(j) = y(k)
      k = k + 1
      pfc(j) = y(k)
      k = k + 1
      hgapt(j) = y(k)
      do i=1,nf+nc
         k = k + 1
         tem(i,j) = y(k)
      end do
      if(inomech .eq. 0)then
         do i=1,nf
            k = k + 1
            efs(i,j) = y(k)
         end do
         do i=1,nf
            k = k + 1
            eft(i,j) = y(k)
         end do
         if(ifcreep .eq. 1)then
            do i=1,nf
               k = k + 1
               efce(i,j) = y(k)
            end do
            do i=1,nf
               k = k + 1
               efch(i,j) = y(k)
            end do
            do i=1,nf
               k = k + 1
               efcr(i,j) = y(k)
            end do
            do i=1,nf
               k = k + 1
               efcz(i,j) = y(k)
            end do
         end if
         do i=1,nf
            k = k + 1
            efh(i,j) = y(k)
         end do
         do i=1,nf
            k = k + 1
            efr(i,j) = y(k)
         end do
         k = k + 1
         efz(j) = y(k)
         do i=1,nf
            k = k + 1
            sigfh(i,j) = y(k)
         end do
         do i=1,nf
            k = k + 1
            sigfr(i,j) = y(k)
         end do
         do i=1,nf
            k = k + 1
            sigfz(i,j) = y(k)
         end do
         do i=1,nc
            k = k + 1
            et(i,j) = y(k)
         end do
         if(iccreep .eq. 1)then
            do i=1,nc
               k = k + 1
               ece(i,j) = y(k)
            end do
            do i=1,nc
               k = k + 1
               ech(i,j) = y(k)
            end do
            do i=1,nc
               k = k + 1
               ecr(i,j) = y(k)
            end do
            do i=1,nc
               k = k + 1
               ecz(i,j) = y(k)
            end do
         end if
         do i=1,nc
            k = k + 1
            eh(i,j) = y(k)
         end do
         do i=1,nc
            k = k + 1
            er(i,j) = y(k)
         end do
         k = k + 1
         ez(j) = y(k)
         do i=1,nc
            k = k + 1
            sigh(i,j) = y(k)
         end do
         do i=1,nc
            k = k + 1
            sigr(i,j) = y(k)
         end do
         do i=1,nc
            k = k + 1
            sigz(i,j) = y(k)
         end do
      end if
   end do
   return
   end
 
!--------------------------------------------------------------------------------------------------
! Read time derivatives from vector y_t
!--------------------------------------------------------------------------------------------------
   subroutine read_from_y_t(y_t)

   integer(c_int) i,j,k
   real(c_double) y_t(maxeq)

   k = 0
   k = k + 1
   dfggen = y_t(k)
   k = k + 1
   dfgrel = y_t(k)
   k = k + 1
!   dgpres = y(k)
   do j=1,nzz
      k = k + 1
      !dqqv1(j) = y_t(k)
      k = k + 1
      dbup(j) = y_t(k)
      k = k + 1
      !dgap(j) = y_t(k)
      k = k + 1
      dpfc(j) = y_t(k)
      k = k + 1
      !dhgapt(j) = y_t(k)
      do i=1,nf+nc
         k = k + 1
         dtem(i,j) = y_t(k)
      end do
      if(inomech .eq. 0)then
         do i=1,nf
            k = k + 1
            defs(i,j) = y_t(k)
         end do
         do i=1,nf
            k = k + 1
!            deft(i,j) = y_t(k)
         end do
         if(ifcreep .eq. 1)then
            do i=1,nf
               k = k + 1
               defce(i,j) = y_t(k)
            end do
            do i=1,nf
               k = k + 1
               defch(i,j) = y_t(k)
            end do
            do i=1,nf
               k = k + 1
               defcr(i,j) = y_t(k)
            end do
            do i=1,nf
               k = k + 1
               defcz(i,j) = y_t(k)
            end do
         end if
         do i=1,nf
            k = k + 1
            defh(i,j) = y_t(k)
         end do
         do i=1,nf
            k = k + 1
!            defr(i,j) = y_t(k)
         end do
         k = k + 1
         defz(j) = y_t(k)
         do i=1,nf
            k = k + 1
!            dsigfh(i,j) = y_t(k)
         end do
         do i=1,nf
            k = k + 1
!            dsigfr(i,j) = y_t(k)
         end do
         do i=1,nf
            k = k + 1
!            dsigfz(i,j) = y_t(k)
         end do
         do i=1,nc
            k = k + 1
!            det(i,j) = y_t(k)
         end do
         if(iccreep .eq. 1)then
            do i=1,nc
               k = k + 1
               dece(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
               dech(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
               decr(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
               decz(i,j) = y_t(k)
            end do
         end if
         do i=1,nc
            k = k + 1
            deh(i,j) = y_t(k)
         end do
         do i=1,nc
            k = k + 1
!            der(i,j) = y_t(k)
         end do
         k = k + 1
         dez(j) = y_t(k)
         do i=1,nc
            k = k + 1
!            dsigh(i,j) = y_t(k)
         end do
         do i=1,nc
            k = k + 1
!            dsigr(i,j) = y_t(k)
         end do
         do i=1,nc
            k = k + 1
!            dsigz(i,j) = y_t(k)
         end do
      end if
   end do
   return
   end

!--------------------------------------------------------------------------------------------------
! Specify absolute error for each variable
!--------------------------------------------------------------------------------------------------
   subroutine abs_err(atol)

   integer(c_int) i,j,k
   real(c_double) atol(maxeq),btol,etol,ftol,gtol,htol,ptol,stol,ttol

!  absolute tolerances:
   btol = 1.0d-4   !burnup (MWd/kg)
   etol = 1.0d-6   !strain (-)
   ftol = 1.0d-7   !fission gas (mol)
   htol = 1.0d-4   !gap conductance (W/m2K)
   gtol = 1.0d-6   !gap (m)
   stol = 1.0d-4   !stress (MPa)
   ttol = 1.0d-6   !temperature (K)

   k = 0
   !fggen
   k = k + 1
   atol(k)=ftol
   !fgrel
   k = k + 1
   atol(k)=ftol
   !gpres
   k = k + 1
   atol(k)=stol
   do j=1,nzz
      !qv(j)
      k = k + 1
      atol(k)=btol
      !bup(j)
      k = k + 1
      atol(k)=btol
      !gap(j)
      k = k + 1
      !pfc(j)
      k = k + 1
      atol(k)=stol
      !hgapt(j)
      k = k + 1
      atol(k)=htol
      do i=1,nf+nc
         !tem
         k = k + 1
         atol(k)=ttol
      end do
      do i=1,nf
         !efs
         k = k + 1
         atol(k)=etol
      end do
      do i=1,nf
         !eft
         k = k + 1
         atol(k)=etol
      end do
      if(ifcreep .eq. 1)then
         do i=1,nf
            !efce
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nf
            !efch
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nf
            !efcr
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nf
            !efcz
            k = k + 1
            atol(k)=etol
         end do
      end if
      do i=1,nf
         !efh
         k = k + 1
         atol(k)=etol
      end do
      do i=1,nf
         !efr
         k = k + 1
         atol(k)=etol
      end do
      !efz
      k = k + 1
      atol(k)=etol
      do i=1,nf
         !sigfh
         k = k + 1
         atol(k)=stol
      end do
      do i=1,nf
         !sigfr
         k = k + 1
         atol(k)=stol
      end do
      do i=1,nf
         !sigfz
         k = k + 1
         atol(k)=stol
      end do
      do i=1,nc
         !et
         k = k + 1
         atol(k)=etol
      end do
      if(iccreep .eq. 1)then
         do i=1,nc
            !ece
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nc
            !ech
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nc
            !ecr
            k = k + 1
            atol(k)=etol
         end do
         do i=1,nc
            !ecz
            k = k + 1
            atol(k)=etol
         end do
      end if
      do i=1,nc
         !eh
         k = k + 1
         atol(k)=etol
      end do
      do i=1,nc
         !er
         k = k + 1
         atol(k)=etol
      end do
      !ez
      k = k + 1
      atol(k)=etol
      do i=1,nc
         !sigh
         k = k + 1
         atol(k)=stol
      end do
      do i=1,nc
         !sigr
         k = k + 1
         atol(k)=stol
      end do
      do i=1,nc
         !sigz
         k = k + 1
         atol(k)=stol
      end do
   end do
   return
   end

end module FRED