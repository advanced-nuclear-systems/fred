!==================================================================================================
! Main driver
!==================================================================================================
program main

   use, intrinsic :: iso_c_binding
   use globals
   use fsundials_core_mod            ! Access SUNDIALS core types, data structures, etc.
   use fida_mod                      ! Fortran interface to IDA
   use fnvector_serial_mod           ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
   use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
   use DAE                           ! ODE functions

   implicit none
   real(c_double) :: t, tout, tret(1), h(1)
   integer(c_int) :: iout, retval, retvalr, nrtfn, i, rootsfound(maxz), k

   integer*8 maxsteps
   integer*4 nnn,maxnef,maxcor

   type(N_Vector), pointer :: sunvec_y                ! sundials solution vector
   type(N_Vector), pointer :: sunvec_yp               ! sundials derivative vector
   type(N_Vector), pointer :: sunvec_av               ! sundials tolerance vector
   type(SUNMatrix), pointer :: sunmat_A               ! sundials matrix
   type(SUNLinearSolver), pointer :: sunlinsol_LS     ! sundials linear solver
   type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
   type(c_ptr) :: ida_mem                             ! IDA memory
   type(c_ptr) :: sunctx                              ! SUNDIALS simulation context

   real(c_double) :: yval(maxeq), ypval(maxeq), rtol, avtol(maxeq)

   retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

!  read input file
   call rdinp()

!  initialize vector of unknowns y
   call write_to_y(yval)

!  initialize vector of time derivatives
   do i = 1,neq
      ypval = 0.0d0
   end do

!  set relative tolerance
   rtol = 1.0d-8
!  set vector of absolute tolerances avtol
   call abs_err(avtol)

!  create serial vectors
   sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
   if (.not. associated(sunvec_y)) then
      print *, 'ERROR: sunvec = NULL'
      stop
   end if

   sunvec_yp => FN_VMake_Serial(neq, ypval, sunctx)
   if (.not. associated(sunvec_yp)) then
      print *, 'ERROR: sunvec = NULL'
      stop
   end if

   sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
   if (.not. associated(sunvec_av)) then
      print *, 'ERROR: sunvec = NULL'
      stop
   end if

!  Call FIDACreate and FIDAInit to initialize IDA memory
   ida_mem = FIDACreate(sunctx)
   if (.not. c_associated(ida_mem)) then
      print *, 'ERROR: ida_mem = NULL'
      stop
   end if

!  Provides required problem and solution specifications, allocates internal memory, and initializes idas
   retval = FIDAInit(ida_mem, c_funloc(residuals), t, sunvec_y, sunvec_yp)
   if (retval /= 0) then
      print *, 'Error in FIDAInit, retval = ', retval
      stop
   end if

!  Call FIDASVtolerances to set tolerances
   retval = FIDASVtolerances(ida_mem, rtol, sunvec_av)
   if (retval /= 0) then
      print *, 'Error in FIDASVtolerances, retval = ', retval
      stop
   end if

!  Create dense SUNMatrix for use in linear solvers
   sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
   if (.not. associated(sunmat_A)) then
      print *, 'ERROR: sunmat = NULL'
      stop
   end if

!  Create dense SUNLinearSolver object
   sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
   if (.not. associated(sunlinsol_LS)) then
      print *, 'ERROR: sunlinsol = NULL'
      stop
   end if

!  Attach the matrix and linear solver
   retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
   if (retval /= 0) then
      print *, 'Error in FIDASetLinearSolver, retval = ', retval
      stop
   end if

!  Create Newton SUNNonlinearSolver object. IDA uses a Newton SUNNonlinearSolver by default, so it is not necessary
!  to create it and attach it. It is done here solely for demonstration purposes.
   sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
   if (.not. associated(sunnonlin_NLS)) then
      print *, 'ERROR: sunnonlinsol = NULL'
      stop
   end if

!  Attach the nonlinear solver
   retval = FIDASetNonlinearSolver(ida_mem, sunnonlin_NLS)
   if (retval /= 0) then
      print *, 'Error in FIDASetNonlinearSolver, retval = ', retval
      stop
   end if

!  Set initial time step
   retval = FIDASetInitStep(ida_mem, 1.0d0)
   if (retval /= 0) then
      print *, 'Error in FIDASetInitStep, retval = ', retval
      stop
   end if

!  Call FIDARootInit to specify the root function with nzz components
   retval = FIDARootInit(ida_mem, nzz, c_funloc(root))
   if (retval /= 0) then
      print *, 'Error in FIDARootInit, retval = ', retval, '; halting'
      stop
   end if

!  save initial outfrd
   call outfrd(0.0d0,0)

   maxsteps = 1000  ! Default (500)
   retval = FIDASetMaxNumSteps(ida_mem, maxsteps)
   maxnef = 1000 ! Default (10)
   retval = FIDASetMaxErrTestFails(ida_mem, maxnef)
   maxcor = 1000 ! Default (4)
   retval = FIDASetMaxNonlinIters(ida_mem, maxcor)
!
!  set maximum timestep
   retval = FIDASetMaxStep(ida_mem, hmax)

!  set initial time
   t = 0.d0

   iout = 0
   nnn = 0
   do while(t .le.tend)

     t = t + dtout

     retval = FIDASolve(ida_mem, t, tret, sunvec_y, sunvec_yp, IDA_NORMAL)

     if (retval .eq. IDA_ROOT_RETURN) then
        retval = FIDAGetRootInfo(ida_mem, rootsfound)
        if (retval < 0) then
           print *, 'Error in FIDAGetRootInfo, retval = ', retval
           stop
        end if
        k = 0
        do i=1,nzz
           k = k + 1
           if(rootsfound(k) .eq. -1)then
              write(*,*)'Gap closed at axial layer ', i
           end if
        end do
     end if
     if (retval < 0) then
        print *, 'Error in FIDASolve, retval = ', retval
        stop
     end if
     nnn = nnn + 1
     call outfrd(t,nnn)

     !retval = FIDAGetCurrentStep(ida_mem, h)
     !dtout = h(1)
  end do

! free memory
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)

end program main


!==================================================================================================
! Calculation of the clad specific heat J/(kg*K)
!     input arguments :
!        tk .... temperature (k)
!==================================================================================================
      real*8 function ccp (tk,cmat)
      implicit none
      real*8 tk,temp(16),cp(16),temp2(17),cp2(17),tf
      character*6 cmat

      if(cmat.eq.'aim1')then
!        L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
         ccp = 431.0d0 + 0.177d0*tk + 8.72d-5*tk**2
      else
         write(*,*)'ccp: wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Returns clad material effective creep rate (1/s)
!     input arguments :
!        tk -- temperature (K)
!        nflux -- neutron flux (n/cm^2*s)
!        mnen -- mean neutron energy (MeV)
!        sg -- effective stress (MPa)
!        cmat -- clad material
!==================================================================================================
      real*8 function ccreep(tk,sg,cmat)
      implicit none
      real*8 tk,sg,tt,E,nflux
      character*6 cmat
      if(cmat.eq.'aim1')then
!        L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
         tt = 1.986d0*tk
!        thermal creep (%/h)
         ccreep = 2.3d14*dexp(-8.46d4/tt)*dsinh(39.72d0*sg/tt)
!        mean neutron energy (MeV)
         E = 2.0d0
!        neutron flux (n/cm2-s)
         nflux = 1.0d18
!        irradiation creep (%/h)
         ccreep = ccreep + 3.2d-24*E*nflux*sg
!        convert from %/h to 1/s
         ccreep=ccreep/100.d0/3.6d3
      else
         ccreep=0.0d0
      end if
      return
      end

!==================================================================================================
! Calculation of the clad Young's modulus
!                unit: MPa
!     input arguments :
!        tk .... temperature (K)
!==================================================================================================
      real*8 function celmod(tk,cmat)
      implicit none
      real*8 tk,tc
      character*6 cmat

      celmod = 0.0d0
      if(cmat.eq.'aim1')then
         tc=tk-273.15d0
         celmod=2.027d11-8.167d7*tc
      else
         write(*,*)'celmod: wrong clad material:',cmat
         stop
      end if
      celmod = celmod/1.0d6
      return
      end

!==================================================================================================
! Calculation of the clad thermal conductivity
!                unit: w/(m*k)
!     input arguments :
!        tk .... temperature (K)
!==================================================================================================
      real*8 function clamb(tk,cmat)
      implicit none
      real*8 tk, tc
      character*6 cmat

      if(cmat.eq.'aim1')then
         tc=tk-273.15d0
         clamb=13.95d0+1.163d-2*tc
      else
         write(*,*)'clamb: wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of the clad Poisson ratio
!                unitless
!     input arguments :
!        tk .... temperature (k)
!==================================================================================================
      real*8 function cpoir(cmat)
      implicit none
      character*6 cmat

      if(cmat.eq.'aim1')then
         cpoir=0.289d0
      else
         write(*,*)'cpoir: wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculates the clad thermal expansion strain
!==================================================================================================
      real*8 function ctexp(tk,cmat)
      implicit none
      real*8 tk, tc
      character*6 cmat

      if(cmat.eq.'aim1')then
!        L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
         tc = tk - 273.15d0
         ctexp=-3.101d-4 + 1.545d-5*tc + 2.75d-9*tc**2
      else
         write(*,*)'wrong cmat option ',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of fuel specific heat
!                  unit : j/(kg*k)
!       input arguments :
!          tk .... temperature (k)
!==================================================================================================
      real*8 function fcp (tk,pucont,fmat)
      implicit none
      real*8 tk,pucont,polp,ucp_,pucp_,tPopov(29),cpuo2Popov(29),cppuo2Popov(29),tliq,tf
      character*3 fmat
!     S.G. Popov, J.J. Carbajo, V.K. Ivanov, G.L. Yoder. Thermophysical properties of MOX and UO2 fuels including the effects of irradiation. ORNL/TM-2000/351 http://web.ornl.gov/~webworks/cpr/v823/rpt/109264.pdf
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
!        melting point
         tliq=3120.0d0-388.1d0*pucont-30.4d0*pucont**2
         if(tk.ge.tliq-1.0d0)then
            fcp=fcp+1.d5*abs(tk-tliq+1.0d0)
         end if
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Returns fuel material effective creep rate (1/s)
!==================================================================================================
      real*8 function fcreep(i,j)
      use globals

      implicit none
      integer i,j
      real*8 polp
      real*8 rate,sg

      fcreep=0.0d0
      if(fmat .eq. 'mox' .and. sigf(i,j) .gt. 0.0d0)then
         sg = sigf(i,j)*1.0d6
!        fission rate per unit volume (fiss/m3-s): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
         rate = qqv1(j)/(1.60214d-13*210.0d0)
         fcreep = fcreepc(1) * sg**fcreepc(2) * rate * dexp(-fcreepc(3)/(rmu*tem(i,j)))
      end if
      return
      end

!==================================================================================================
! Calculation of the fuel Young's modulus, unit: MPa
!     input arguments :
!        tk .... temperature (k)
!==================================================================================================
      real*8 function felmod(tk,den,cont,fmat)
      implicit none
      real*8 tk,den,cont,tden,por
      character*6 fmat

      if(fmat.eq.'mox')then
         tden=11460.d0*cont+10960.d0*(1.d0-cont)
         por=1.d0-den/tden
!        MATPRO correlation
         felmod=(2.334d11*(1.0d0-2.752d0*por)*(1.d0-1.0915d-4*tk))*(1.d0+0.15d0*cont)
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      felmod = felmod/1.0d6
      return
      end

!==================================================================================================
! Calculation of fission gas release
!==================================================================================================
      real*8 function fgr(j)

      use globals
      implicit none
      integer j,i
      real*8 a,fgrR,fgrU,aR,r

!     Waltar and Reynolds model
      if(bup(j) .eq. 0.0d0)then
         fgr = 0.0d0
      else
!        calculate release fraction for the restructured fuel
         a = 4.7d0/bup(j) * (1.0d0-dexp(-bup(j)/5.9d0))
         fgrR = max(0.0d0, 1.0d0 - a)
      
!        calculate release fraction for the unrestructured fuel
         if(bup(j) .le. 3.5d0)then
             fgrU = 0.0d0
         else    
             a = 25.6d0/(bup(j)-3.5d0) * (1.0d0-dexp(-(bup(j)/3.5d0-1.0d0)))*dexp(-0.0125d0*(ql(j)/1.d3))
             if(bup(j) .ge. 49.2d0) a = a*exp(-0.3d0*(bup(j)-49.2d0))
             fgrU = max(0.0, 1.0d0 - a)
         end if
!        fractional fuel areas associated with restructured zone
         aR = 0.0d0
!        assume that the temperature 1500 C defines the boundary between unrestructured and restructured zones (Waltar & Reynolds, p.198).
         do i=1,nf
            if(tem(i,j)-273.15d0 .ge. 1500.0d0)aR = aR + az0(i)
         end do
         aR = aR/azf0
         fgr = (1.0d0 - aR)*fgrU + aR*fgrR
      end if
      
      return
      end

!==================================================================================================
! Calculation of the fuel local thermal conductivity, unit : w/(m*k)
!     input arguments :
!        dens - density (kg/m3)
!        bup - burn-up (mwd/kgu)
!        tk - temperature (k)
!        pucon - Pu fraction (-)
!        sto - fuel stoichiometry
!        fmat - fuel material index
!==================================================================================================
      real*8 function flamb(dens,bup,tk,pucon,sto,fmat)
      implicit none
      real*8 dens,bup,tk,pucon,sto,ac,tden,por
      character*3 fmat

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
      end

!==================================================================================================
! Calculation of the fuel Poisson ratio, unitless
!==================================================================================================
      real*8 function fpoir(fmat)
      implicit none
      character*6 fmat

!     MATPRO model
      if(fmat.eq.'mox')then
         fpoir=0.276d0
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Volumetric fuel swelling rate (-)
!==================================================================================================
      real*8 function fswel(i,j)
      use globals
      implicit none
      integer i,j
      real*8 rate, dT, E, B

      fswel=0.0d0
      if(ifswel .eq. 1)then
         if(fmat .eq. 'mox')then
!           !fission rate per unit volume (fiss/m3-s): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
            rate = qqv1(j)/(1.60214d-13*210.0d0)
!           MATPRO model for swelling due to gaseous FPs
            if(tem(i,j) .lt. 2800.0d0)then
               dT = 2800.0d0-tem(i,j)
!              total energy generated per unit volume (J/m3)
               E = bup(j)*1.0d6*8.64d4*rof0
!              total number of fissions per unit volume (fiss/m3): 1 J/fiss = 210 MeV/fiss * 1.60214d-13 J/MeV
               B = E/(1.60214d-13*210.0d0)
               fswel = 8.8d-56 * dT**11.73d0 * exp(-0.0162d0*dT) * exp(-8.0d-27*B) * rate
            else
               fswel = 0.0d0
            end if

!           MATPRO model for swelling due to solid FPs
            fswel = fswel + 2.5d-29*rate
         end if
      end if
      return
      end

!==================================================================================================
! Calculation of fuel thermal expansion, unit : m/m
!       input arguments :
!          tk .... temperature (k)
!==================================================================================================
      real*8 function ftexp(tk,cont,fmat)
      implicit none
      real*8 tk,cont,tc,ftexpUO2,ftexpPuO2
      character*3 fmat

      if(fmat.eq.'mox')then
!        MATPRO model
         ftexpUO2  = 1.0d-5*tk - 3.0d-3 + 4.d-2*dexp(-5000.d0/tk)
         ftexpPuO2 = 0.9d-5*tk - 2.7d-3 + 7.d-2*dexp(-5072.d0/tk)
         ftexp = ftexpUO2*(1.0-cont) + ftexpPuO2*cont 
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Open and closed fuel/clad gap conductance
!  j - axial slice
!==================================================================================================
      real*8 function gaphtc(j)

      use globals
      implicit none
      integer j,i,jj
      real*8 gtemp,fe,sbc,emissf,emissc,cmhard,r,conf,conc,fkm,flamb,clamb,ahe,pgas,gask_(4), &
     &       xmol_(4),sumx,tf,tc,hgap(3),mw(4),gap_
!             he        ar        kr      xe
      data mw/4.0026d0, 39.948d0, 83.8d0, 131.30d0/

!     Stefan-Boltzmann constant (w/m**2-k**4)
      data sbc/5.6697d-8/

      tf=tem(nf,j)
      tc=tem(nf+1,j)

!     average gap temperature
      gtemp=(tf+tc)/2.d0
!     helium thermal conductivity
      gask_(1)=2.639d-3*gtemp**0.7085d0
!     argon thermal conductivity
      gask_(2)=2.986d-4*gtemp**0.7224d0
!     kripton thermal conductivity
      gask_(3)=8.247d-5*gtemp**0.8363d0
!     xenon thermal conductivity
      gask_(4)=4.351d-5*gtemp**0.8618d0
      if(gmat.eq.'he')then
         xmol_(1)=mu0
         xmol_(2)=0.d0
      else ! if(gmat.eq.'ar')then
         xmol_(1)=0.d0
         xmol_(2)=mu0
      end if

!     number of moles of released fission gases (0.8846xe,0.0769kr,0.0385he)
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

!     open gap conductance
      ahe=0.425d0-2.3d-4*dmin1(gtemp,1000.d0)
      pgas=dmax1(gpres,0.1d0)*1.0d6
!     accomodation distance
      ajump(j)=0.024688d0*gask(j)*dsqrt(gtemp)*2.d0/(pgas*ahe)

!     using FRAPCON model reduce the gap width due to fuel cracking while calculating the gap conductance 
      if(ifreloc .eq. 1)then
         reloc(j)=0.25d0 + dmax1( 5.0d0, dmin1(10.d0,ql(j)/4000.d0) )*( dmin1( 1.0d0, bup(j)/5.d0) + 1.0d0 )/100.d0
      else
         reloc(j)=0.0d0
      end if

!     effective roughness
      r=sqrt(rufc**2+ruff**2)
      gap_=max(gap(j), r)
 
      gapth(j) = gap_*(1.0d0 - reloc(j)) + ajump(j)
      hgap(1)=gask(j)/gapth(j)

!     radiation conductance
!     set fuel and clad surface emissivities to reasonable value.
!     according to KIT recommendation (Struwe+Schikorr) 15.06.2011 
      emissf = 0.8d0
      emissc = 0.9d0

      fe=1.0d0/(1.0d0/emissf+(rfo0/rci0)*(1.0d0/emissc-1.0d0))
      hgap(2)=sbc*fe*(tf**2+tc**2)*(tf+tc)

      if(flag(j).eq.'open')then
         hgap(3)=0.d0
      else
!        average cladding temperature
         tc = 0.0d0
         do i=1,nc
            tc = tc + tem(nf+i,j)*az0(nf+i)
         end do
         tc = tc/azc0
!        Meyer hardness for SS (from SAS4A hgap routine)
         if(tc.le.893.9203d0)then
            cmhard=5.961d9*tc**(-0.206d0)
         else
            cmhard=2.750d28*tc**(-6.530d0)
         end if
         if(cmhard.lt.1.d5)cmhard=1.d5

!        effective roughness
         r=sqrt(rufc**2+ruff**2)

!        effective conductivity
         conf=flamb(rof0,bup(j),tf,pucont,sto0,fmat)
         conc=clamb(tc,cmat)
         fkm=2.d0*conf*conc/(conf+conc)

!        closed gap conductance
         hgap(3)=33.3d0*fkm*(pfc(j)-gpres)*1.0d6/(cmhard*sqrt(r))
      end if
      gaphtc = hgap(1) + hgap(2) + hgap(3)

      hgapi(1,j)=hgap(1)
      hgapi(2,j)=hgap(2)
      hgapi(3,j)=hgap(3)
      
      return
      end

!==================================================================================================
! Internal gas pressure calculation
!==================================================================================================
      real*8 function gpresf(time)

      use globals
      implicit none
      integer j
      real*8 time,v_t,tgap,thol,tp,polp

!     table of gas plenum temperature (K) vs time
      if(ntple.gt.0)then
         tp = polp(time,tple,ttple,ntple)
      else
         tp = tem0
      end if

!     volume-to-temperature ratio
      v_t = vgp/tp
      do j=1,nzz
!       gap and hole temperature
        tgap = (tem(nf,j) + tem(nf+1,j))/2.d0
        thol = tem(1,j)
!       gap vol/tem + hole vol/tem
        v_t = v_t + pi*dz0(j)*((rfo0**2 - rci0**2)/tgap + rfi0**2/thol)
      end do

!     free volume gas pressure (MPa)
      gpresf = (mu0 + fgrel)*rmu / v_t * 1d-6

      return
      end

!==================================================================================================
! Write outfrd file
!==================================================================================================
      subroutine outfrd(time,nnn)

      use globals
      implicit none
      integer nnn
      integer i,j
      real*8 time
      real*8 flamb
      character*12 s
      character*20 title

      write(s,'(i12.12)') nnn
      open(700,file='outfrd' // s)
      title = "time (s)            "
      write(700,2002)title,time
      
      title = "time (d)            "
      write(700,2002)title,time/86400.d0
      
      write(700,'(500a)')('-',j=1,100)
      title = "rof0 (kg/m3)        "
      write(700,2001)title,rof0
      
      title = "fggen (cm3)         "
      write(700,2001)title,fggen*rmu*293.15d0/1.0d5 * 1.0d6
      
      title = "fgrel (cm3)         "
      write(700,2001)title,fgrel*rmu*293.15d0/1.0d5 * 1.0d6
      
      title = "fgrel (%)           "
      write(700,2001)title,fgrel/(fggen+1.d-21)*100.d0
      
      title = "gpres (MPa)         "
      write(700,2001)title,gpres
      
      write(700,'(500a)')('-',j=1,100)
      title = "z (m)               "
      write(700,2001)title,(dz0(j)*(real(j)-0.5d0),j=1,nzz)
      
      title = "tfin (C)            "
      write(700,2001)title,(tem(1,j)-273.15d0,j=1,nzz)
      
      title = "tfout (C)           "
      write(700,2001)title,(tem(nf,j)-273.15d0,j=1,nzz)
      
      title = "tcin (C)            "
      write(700,2001)title,(tem(nf+1,j)-273.15d0,j=1,nzz)
      
      title = "tcout (C)           "
      write(700,2001)title,(tem(nf+nc,j)-273.15d0,j=1,nzz)
      
      title = "rfi (m)             "
      write(700,2001)title,(rfi(j),j=1,nzz)
      
      title = "rfo (m)             "
      write(700,2001)title,(rfo(j),j=1,nzz)
      
      title = "rci (m)             "
      write(700,2001)title,(rci(j),j=1,nzz)
      
      title = "rco (m)             "
      write(700,2001)title,(rco(j),j=1,nzz)
      
      title = "dzf (m)             "
      write(700,2001)title,(dzf(j),j=1,nzz)
      
      title = "dzc (m)             "
      write(700,2001)title,(dzc(j),j=1,nzz)
      
      title = "reloc (-)           "
      write(700,2001)title,(reloc(j),j=1,nzz)
      
      title = "qv (W/m3)           "
      write(700,2001)title,(qqv1(j),j=1,nzz)
      
      title = "ql (W/cm)           "
      write(700,2001)title,(ql(j)/100.d0,j=1,nzz)
      
      title = "ql2 (W/cm)          "
      write(700,2001)title,(ql(j)/100.d0*dz0(j)/dzf(j),j=1,nzz)
      
      title = "bup (MWd/kg)        "
      write(700,2001)title,(bup(j),j=1,nzz)
      
      title = "gap (m)             "
      write(700,2001)title,(gap(j),j=1,nzz)
      
      title = "gapth (m)             "
      write(700,2001)title,(gapth(j),j=1,nzz)
      
      title = "hgap (W/m2K)        "
      write(700,2001)title,(hgapt(j),j=1,nzz)
      
      title = "hgap1 (W/m2K)       "
      write(700,2001)title,(hgapi(1,j),j=1,nzz)
      
      title = "hgap2 (W/m2K)       "
      write(700,2001)title,(hgapi(2,j),j=1,nzz)
      
      title = "hgap3 (W/m2K)       "
      write(700,2001)title,(hgapi(3,j),j=1,nzz)
      
      title = "ajump (m)           "
      write(700,2001)title,(ajump(j),j=1,nzz)
      
      title = "gask (W/m-K)        "
      write(700,2001)title,(gask(j),j=1,nzz)
      
      title = "pfc (MPa)           "
      write(700,2001)title,(pfc(j),j=1,nzz)
      
      write(700,'(500a)')('-',j=1,100)
      title = "r (m)               "
      write(700,2001)title,(rad0(i),i=1,nf+nc)
      
      title = "temperature (C)     "
      do j=1,nzz
         write(700,2000)title,j,(tem(i,j)-273.15d0,i=1,nf+nc)
      end do
      
      title = "sig h (MPa)         "
      do j=1,nzz
         write(700,2000)title,j,(sigfh(i,j),i=1,nf),(sigh(i,j),i=1,nc)
      end do
      
      title = "sig r (MPa)         "
      do j=1,nzz
         write(700,2000)title,j,(sigfr(i,j),i=1,nf),(sigr(i,j),i=1,nc)
      end do
      
      title = "sig z (MPa)         "
      do j=1,nzz
         write(700,2000)title,j,(sigfz(i,j),i=1,nf),(sigz(i,j),i=1,nc)
      end do
      
      title = "eps h total (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efh(i,j)*100.d0,i=1,nf),(eh(i,j)*100.d0,i=1,nc)
      end do
      
      title = "eps r total (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efr(i,j)*100.d0,i=1,nf),(er(i,j)*100.d0,i=1,nc)
      end do
      
      title = "eps z total (%)    "
      do j=1,nzz
         write(700,2000)title,j,efz(j)*100.d0,ez(j)*100.d0
      end do

      title = "eps r creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efcr(i,j)*100.0d0,i=1,nf),(ecr(i,j)*100.0d0,i=1,nc)
      end do
      
      title = "eps h creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efch(i,j)*100.0d0,i=1,nf),(ech(i,j)*100.0d0,i=1,nc)
      end do
      
      title = "eps z creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efcz(i,j)*100.0d0,i=1,nf),(ecz(i,j)*100.0d0,i=1,nc)
      end do
      
      title = "eps thermal lin (%) "
      do j=1,nzz
         write(700,2000)title,j,(eft(i,j)*100.0d0,i=1,nf),(et(i,j)*100.0d0,i=1,nc)
      end do
      
      title = "eps swell vol (%)   "
      do j=1,nzz
         write(700,2000)title,j,(efs(i,j)*100.0d0,i=1,nf)
      end do

      close(700)
2000  format(a20,1x,"iz:",i2,500(1x,1pe12.5))
2001  format(a20,1x,"   ",2x,500(1x,1pe12.5))
2002  format(a20,1x,"   ",2x,500(1x,1pe12.5))

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

!==================================================================================================
! polp is a linear interpolation function whose value is equal to
! interpolated value.
! tc-is the value of the independent variable
! crf-dependent variable array
! crft-independent variable array
! n-number of points in variable arrays
!==================================================================================================
      function polp (tc,crf,crft,n)

      implicit real*8 (a-h,o-z)
      dimension crf(n),crft(n)
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
      end

!==================================================================================================
! Reads input deck and assigns initial values
!==================================================================================================
      subroutine rdinp()

      use globals
      implicit none
      integer numlin,i1,i2,i,j
      real*8 ftexp,felmod,fpoir,gpresf,ctexp,celmod,cpoir,gaphtc
      real*8 r1,eft0,fyng0,fpoir0,ef0,sigf0,et0,cyng0,cpoir0,tfree0,vfree0
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

               else if(w2 .eq. 'CLAD_CREEP')then
                  iccreep=1

               else if(w2 .eq. 'FUEL_RELOC')then
                  ifreloc=1

               else if(w2 .eq. 'NOMECH')then
                  inomech=1

               end if   

!           initial temperature card
            else if(i1.eq.000002) then
               backspace 700
               read(700,*,err=1000)i1,tem0

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
!              cmat - clad material (aim1)
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
      eft0 = ftexp(tem0,pucont,fmat) - ftexp(293.15d0,pucont,fmat)
      fyng0 = felmod(tem0,rof0,pucont,fmat)
      fpoir0 = fpoir(fmat)
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
      et0 = ctexp(tem0,cmat) - ctexp(293.15d0,cmat)
      cyng0 = celmod(tem0,cmat)
      cpoir0 = cpoir(cmat)
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

!==================================================================================================
! Update geometry for axial layer j
!==================================================================================================
   subroutine update_geom(j)

   use globals
   implicit none
   integer i,j,k
   real*8 efr_(maxr),er_(maxr)

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

!==================================================================================================
! Write variables to vector y
!==================================================================================================
      subroutine write_to_y(y)

      use globals
      implicit none
      integer i,j,k
      real*8, intent(inout) :: y(maxeq)

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

!==================================================================================================
! Read variables from vector y
!==================================================================================================
      subroutine read_from_y(y)

      use globals
      implicit none
      integer i,j,k
      real*8 y(maxeq)

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
 
!==================================================================================================
! Read time derivatives from vector y_t
!==================================================================================================
      subroutine read_from_y_t(y_t)

      use globals
      implicit none
      integer i,j,k
      real*8 y_t(maxeq)

      k = 0
      k = k + 1
      dfggen = y_t(k)
      k = k + 1
      dfgrel = y_t(k)
      k = k + 1
!      dgpres = y(k)
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
!               deft(i,j) = y_t(k)
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
!               defr(i,j) = y_t(k)
            end do
            k = k + 1
            defz(j) = y_t(k)
            do i=1,nf
               k = k + 1
!               dsigfh(i,j) = y_t(k)
            end do
            do i=1,nf
               k = k + 1
!               dsigfr(i,j) = y_t(k)
            end do
            do i=1,nf
               k = k + 1
!               dsigfz(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
!               det(i,j) = y_t(k)
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
!               der(i,j) = y_t(k)
            end do
            k = k + 1
            dez(j) = y_t(k)
            do i=1,nc
               k = k + 1
!               dsigh(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
!               dsigr(i,j) = y_t(k)
            end do
            do i=1,nc
               k = k + 1
!               dsigz(i,j) = y_t(k)
            end do
         end if
      end do
      return
      end

!==================================================================================================
! Specify absolute error for each variable
!==================================================================================================
      subroutine abs_err(atol)

      use globals
      implicit none
      integer i,j,k
      real*8 atol(maxeq),btol,etol,ftol,gtol,htol,ptol,stol,ttol

!     absolute tolerances:
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
