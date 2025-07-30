!==================================================================================================
! Main driver
!==================================================================================================
program main

   use, intrinsic :: iso_c_binding
   use FRED
   use globals
   use fsundials_core_mod            ! Access SUNDIALS core types, data structures, etc.
   use fida_mod                      ! Fortran interface to IDA
   use fnvector_serial_mod           ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
   use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver

   implicit none
   real(c_double) :: t, tout, tret(1), h(1)
   integer(c_int) :: iout, retval, retvalr, nrtfn, i, j, rootsfound(2*maxz+maxtab), k

   integer(c_long) maxsteps
   integer(c_int) nnn,maxnef,maxcor

   type(N_Vector), pointer :: sunvec_y                ! sundials solution vector
   type(N_Vector), pointer :: sunvec_yp               ! sundials derivative vector
   type(N_Vector), pointer :: sunvec_av               ! sundials tolerance vector
   type(SUNMatrix), pointer :: sunmat_A               ! sundials matrix
   type(SUNLinearSolver), pointer :: sunlinsol_LS     ! sundials linear solver
   type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
   type(c_ptr) :: ida_mem                             ! IDA memory
   type(c_ptr) :: sunctx                              ! SUNDIALS simulation context

   real(c_double) :: yval(maxeq), ypval(maxeq), avtol(maxeq)
   real(c_double) :: next_recommended_step(1)

   retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

!  read input file
   call rdinp()

   if(irst .eq. -1)then
!     set initial time
      t = 0.d0
!     initialize vector of unknowns y
      call write_to_y(yval)
!     initialize vector of time derivatives
      do i = 1,neq
         ypval = 0.0d0
      end do
!     time step number
      nnn = 0
   else
      nnn = irst
!     read initial time, restart time step number nnn, vector of unknowns y, and vector of time derivatives yp from rstfrd
      call rstfrd('read ',t,nnn,yval,ypval)
   end if

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
   retval = FIDARootInit(ida_mem, nzz + ntfix, c_funloc(root))
   if (retval /= 0) then
      print *, 'Error in FIDARootInit, retval = ', retval, '; halting'
      stop
   end if

   maxsteps = 50000  ! Default (500)
   retval = FIDASetMaxNumSteps(ida_mem, maxsteps)
   maxnef = 1000 ! Default (10)
   retval = FIDASetMaxErrTestFails(ida_mem, maxnef)
   maxcor = 1000 ! Default (4)
   retval = FIDASetMaxNonlinIters(ida_mem, maxcor)
!
!  set maximum timestep
   retval = FIDASetMaxStep(ida_mem, hmax)

   iout = 0
   if(irst .eq. -1)then
!     save initial outfrd and rstfrd
      call outfrd(t,nnn)
      call rstfrd('write',t,nnn,yval,ypval)
   end if

   do while(t .le.tend)

     t = t + dtout

     retval = FIDASolve(ida_mem, t, tret, sunvec_y, sunvec_yp, IDA_NORMAL)

     if(retval .eq. IDA_ROOT_RETURN) then
        retval = FIDAGetRootInfo(ida_mem, rootsfound)
        if (retval < 0) then
           print *, 'Error in FIDAGetRootInfo, retval = ', retval
           stop
        end if
        k = 0
!       gap closure event
        do j=1,nzz
           k = k + 1
           if(rootsfound(k) .eq. -1)then
              write(*,*)'Gap closed at axial layer ', j
           end if
        end do

     else if (retval < 0) then
        print *, 'Error in FIDASolve, retval = ', retval
        stop
     end if
     nnn = nnn + 1
     call outfrd(t,nnn)
     call rstfrd('write',t,nnn,yval,ypval)
     
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
