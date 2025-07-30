module globals
   use, intrinsic :: iso_c_binding
   implicit none
!=======================================================================
!  ARRAY DIMENSIONS:
   integer, parameter :: maxeq=100000           !maximum number of equations
   integer, parameter :: maxz=30                !maximum number of axial nodes
   integer, parameter :: maxr=30                !maximum number of radial nodes
   integer, parameter :: maxtab=100             !maximum number of entries in tables

!=======================================================================
!  DIMENSIONS OF TABLES:
   integer nqv                                  !number of entries in table of fuel rod power density vs time
   integer nfsto                                !number of fuel stoichiometry vs fuel burnup table entries
   integer ntco                                 !number of fuel rod clad outer temperature table entries
   integer ntple                                !number of gas plenum temperature table entries
   integer ntfix                                !number of fixed times

!=======================================================================
!  FLAGS:
   integer ifgr                                 !option of fission gas release calculation (0-no, 1-yes)
   integer ifcreep                              !option of fuel creep calculation (0-no, 1-yes)
   integer ifswel                               !option of fuel swelling calculation (0-no, 1-yes)
   integer iccreep                              !option of clad creep calculation (0-no, 1-yes)
   integer icplas                               !option of clad plastic calculation (0-no, 1-yes)
   integer ifreloc                              !option of fuel relocation calculation (0-no, 1-yes)
   integer inomech                              !option of no mechanics calculation (0-no, 1-yes)
   integer irst                                 !rstfrd number to read from for restrat

   character*4 gstate(maxz)                     !flag of gap state ('clos', 'open')

!=======================================================================
!  EQUATION PARAMETERS:
   integer*8 neq                                !number of equations
   integer*8 neq_j                              !number of equations for an axial layer
   type ideq                                    !type equation id
      character(len=10):: id                    ! string
      integer :: j                              ! axial index
      integer :: i                              ! radial index
   end type                                     !type end
   type (ideq) :: id_eq(maxeq)

!=======================================================================
!  TIME PARAMETERS:
   real(c_double) tend                          !terminal time (s)
   real(c_double) dtout                         !output time step (s)
   real(c_double) hmax                          !maximum time step (s)
   real(c_double) tfix(maxtab)                  !fixed times (s)

!=======================================================================
!  TOLERANCE PARAMETERS:
   real(c_double) rtol                          !relative tolerance
   real(c_double) btol                          !absolute tolerance for burnup (MWd/kgHM)
   real(c_double) etol                          !absolute tolerance for strain (m/m)
   real(c_double) ftol                          !absolute tolerance for fission gas amount (mol)
   real(c_double) htol                          !absolute tolerance for gap conductance (W/m2K)
   real(c_double) gtol                          !absolute tolerance for gap width (m)
   real(c_double) stol                          !absolute tolerance for stress (MPa)
   real(c_double) ttol                          !absolute tolerance for temperature (K)

!=======================================================================
!  GEOMETRY PARAMETERS:
   integer nf                                   !number of radial fuel nodes
   integer nc                                   !number of radial clad nodes
   integer nzz                                  !number of axial layers

   real(c_double) gap(maxz)                     !radial fuel-clad gap (m)
   real(c_double) gapth(maxz)                   !"thermal" gap width (m)
   real(c_double) reloc(maxz)                   !fraction of the gap closed due to fuel relocation
   real(c_double) gask(maxz)                    !gas thermal conductivity (W/mK)
   real(c_double) rad0(maxr)                    !initial node radius (m)
   real(c_double) rad_0(maxr)                   !initial radius of boundary between nodes(m)
   real(c_double) rad(maxr,maxz)                !node radius (m)
   real(c_double) rfi0                          !initial inner fuel radius (m)
   real(c_double) rfi(maxz)                     !inner fuel radius (m)
   real(c_double) rfo0                          !initial outer fuel radius (m)
   real(c_double) rfo(maxz)                     !outer fuel radius (m)
   real(c_double) rci0                          !initial inner clad radius (m)
   real(c_double) rci(maxz)                     !inner clad radius (m)
   real(c_double) rco0                          !initial outer clad radius (m)
   real(c_double) rco(maxz)                     !outer clad radius (m)
   real(c_double) dz0(maxz)                     !initial node height (m)
   real(c_double) dzf(maxz)                     !fuel node height (m)
   real(c_double) dzc(maxz)                     !clad node height (m)
   real(c_double) drf0                          !initial fuel node thickness (m)
   real(c_double) drf(maxr,maxz)                !fuel node thickness (m)
   real(c_double) drc0                          !initial clad node thickness (m)
   real(c_double) drc(maxr,maxz)                !clad node thickness (m)
   real(c_double) az0(maxr)                     !initial node cross sectional area (m2)
   real(c_double) azf0                          !initial fuel cross sectional area (m2)
   real(c_double) azc0                          !initial clad cross sectional area (m2)
   real(c_double) vgp                           !gas plenum volume (m3)
   real(c_double) h0                            !initial fuel height (m)

   character*6 fmat
   character*6 gmat
   character*6 cmat

!=======================================================================
!  TEMPERATURE PARAMETERS:
   real(c_double) tem0                          !initial temperature (K)
   real(c_double) tem(maxr,maxz)                !temperature in node (K)
   real(c_double) dtem(maxr,maxz)               !time derivative of temperature in node (K)
   real(c_double) tco(maxz,maxtab)              !fuel rod clad outer temperature (K) in table of fuel rod clad outer temperature vs time
   real(c_double) ttco(maxtab)                  !time (s) in table of fuel rod clad outer temperature vs time

!=======================================================================
!  ENERGY PARAMETERS:
   real(c_double) qv(maxz,maxtab)               !fuel rod power density (W/m3) in table of fuel rod power density vs time
   real(c_double) tqv(maxtab)                   !time (s) in table of fuel rod power density vs time
   real(c_double) ql(maxz)                      !fuel rod linear power (W/m)
   real(c_double) qqv1(maxz)                    !fuel power density (W/m3)

!=======================================================================
!  GAS PLENUM PARAMETERS:
   real(c_double) gpres0                        !initial inner gas pressure (Pa)
   real(c_double) gpres                         !inner gas pressure (Pa)
   real(c_double) mu0                           !as fabricated gas amount (mol)
   real(c_double) tple(maxtab)                  !gas plenum temperature (K) in table of gas plenum temperature vs time
   real(c_double) ttple(maxtab)                 !time (s) in table of gas plenum temperature vs time

!=======================================================================
!  COOLANT PARAMETERS (FOR BASE IRRADIATION ANALYSIS):
   real(c_double) pcool                         !coolant pressure (Pa)

!=======================================================================
!  MATERIAL PROPERTIES PARAMETERS:
   real(c_double) bup(maxz)                     !fuel burnup (MWd/kgU)
   real(c_double) dbup(maxz)                    !time derivative of fuel burnup (MWd/kgU/s)
   real(c_double) pucont                        !plutonium content in fuel (-)
   real(c_double) sto0                          !initial fuel stoichiometry (-)
   real(c_double) rof0                          !initial fuel density (kg/m3)
   real(c_double) roc0                          !initial clad density (kg/m3)
   real(c_double) rof(maxr,maxz)                !fuel density (kg/m3)
   real(c_double) kfuel(maxr,maxz)              !fuel thermal condcutivity (W/mK)
   real(c_double) kclad(maxr,maxz)              !clad thermal condcutivity (W/mK)
   real(c_double) bsto(maxtab)                  !burnup (MWd/kgHM) in table of fuel stoichiometry vs. burnup
   real(c_double) stob(maxtab)                  !fuel stoichiometry (-) in table of fuel stoichiometry vs. burnup
   real(c_double) fswelmlt                      !fuel swelling rate multiplier
   real(c_double) zrcont                        !Zr content in U-Pu_Zr fuel

!=======================================================================
!  GAP PARAMETERS:
   real(c_double) pfc(maxz)                     !pellet-to-cladding contact pressure (MPa)
   real(c_double) dpfc(maxz)                    !time derivative of pellet-to-cladding contact pressure (MPa/s)
   real(c_double) rufc                          !cladding inner surface roughness (m)
   real(c_double) ruff                          !fuel outer surface roughness (m)
   real(c_double) hgapt(maxz)                   !gap conductance (W/m**2-K)
   real(c_double) hgapi(3,maxz)                 !components of gap conductance (W/m**2-K)
   real(c_double) xmol(5)                       !number of moles of gas components (-)
   real(c_double) ajump(maxz)                   !solid-gas temperature jump (m)

!=======================================================================
!  CONSTANTS:
   real(c_double), parameter :: rmu=8.314d0             !universal gas constant
   real(c_double), parameter :: avo=6.0247d23           !Avogadro number
   real(c_double), parameter :: pi=3.141592653589793d0  !pi number

!=======================================================================
!  MECHANICAL PARAMETERS:
   real(c_double) sigh(maxr,maxz)               !clad hoop stress (Pa)
   real(c_double) dsigh(maxr,maxz)              !clad hoop stress time development (Pa/s)
   real(c_double) sigr(maxr,maxz)               !clad radial stress (Pa)
   real(c_double) sigz(maxr,maxz)               !clad axial stress (Pa)
   real(c_double) sig(maxr,maxz)                !clad effective stress (Pa)
   real(c_double) eh(maxr,maxz)                 !clad total hoop strain (m/m)
   real(c_double) er(maxr,maxz)                 !clad total radial strain (m/m)
   real(c_double) ez(maxz)                      !clad total axial strain (m/m)
   real(c_double) deh(maxr,maxz)                !clad total hoop strain time derivative (1/s)
   real(c_double) dez(maxz)                     !clad total axial strain time derivative (1/s)
   real(c_double) ece(maxr,maxz)                !clad creep effective strain (m/m)
   real(c_double) ech(maxr,maxz)                !clad creep hoop strain (m/m)
   real(c_double) ecr(maxr,maxz)                !clad creep radial strain (m/m)
   real(c_double) ecz(maxr,maxz)                !clad creep axial strain (m/m)
   real(c_double) dece(maxr,maxz)               !time derivative of clad creep effective strain (1/s)
   real(c_double) dech(maxr,maxz)               !time derivative of clad creep hoop strain (1/s)
   real(c_double) decr(maxr,maxz)               !time derivative of clad creep radial strain (1/s)
   real(c_double) decz(maxr,maxz)               !time derivative of clad creep axial strain (1/s)
   real(c_double) epe(maxr,maxz)                !clad plastic effective strain (m/m)
   real(c_double) eph(maxr,maxz)                !clad plastic hoop strain (m/m)
   real(c_double) epr(maxr,maxz)                !clad plastic radial strain (m/m)
   real(c_double) epz(maxr,maxz)                !clad plastic axial strain (m/m)
   real(c_double) depe(maxr,maxz)               !time derivative of clad plastic effective strain (1/s)
   real(c_double) deph(maxr,maxz)               !time derivative of clad plastic hoop strain (1/s)
   real(c_double) depr(maxr,maxz)               !time derivative of clad plastic radial strain (1/s)
   real(c_double) depz(maxr,maxz)               !time derivative of clad plastic axial strain (1/s)
   real(c_double) et(maxr,maxz)                 !clad thermal linear strain (m/m)
   real(c_double) sigfh(maxr,maxz)              !fuel hoop stress (Pa)
   real(c_double) sigfr(maxr,maxz)              !fuel radial stress (Pa)
   real(c_double) sigfz(maxr,maxz)              !fuel axial stress (Pa)
   real(c_double) sigf(maxr,maxz)               !fuel effective stress (Pa)
   real(c_double) efh(maxr,maxz)                !fuel total hoop strain (m/m)
   real(c_double) efr(maxr,maxz)                !fuel total radial strain (m/m)
   real(c_double) efz(maxz)                     !fuel total axial strain (m/m)
   real(c_double) defh(maxr,maxz)               !fuel total hoop strain time derivative (1/s)
   real(c_double) defz(maxz)                    !fuel total axial strain time derivative (1/s)
   real(c_double) eft(maxr,maxz)                !fuel linear thermal expansion (m/m)
   real(c_double) efs(maxr,maxz)                !fuel volumetric swelling strain (m/m)
   real(c_double) defs(maxr,maxz)               !time derivative of fuel volumetric swelling strain (1/s)
   real(c_double) fcreepc(3)                    !user input constants for the fuel creep rate law
   real(c_double) efce(maxr,maxz)               !fuel creep effective strain (m/m)
   real(c_double) efch(maxr,maxz)               !fuel creep hoop strain (m/m)
   real(c_double) efcr(maxr,maxz)               !fuel creep radial strain (m/m)
   real(c_double) efcz(maxr,maxz)               !fuel creep axial strain (m/m)
   real(c_double) defce(maxr,maxz)              !time derivative of fuel creep effective strain (1/s)
   real(c_double) defch(maxr,maxz)              !time derivative of fuel creep hoop strain (1/s)
   real(c_double) defcr(maxr,maxz)              !time derivative of fuel creep radial strain (1/s)
   real(c_double) defcz(maxr,maxz)              !time derivative of fuel creep axial strain (1/s)

!=======================================================================
! FISSION GAS RELEASE AND SWELLING PARAMETERS:
   real(c_double) fggen                         !amount of fission gas produced in fuel rod (mol/s)
   real(c_double) fgrel                         !amount of fission gas released in fuel rod (mol)
   real(c_double) dfggen                        !time derivative of amount of fission gas produced in fuel rod (mol/s)
   real(c_double) dfgrel                        !time derivative of amount of fission gas released in fuel rod (mol/s)

end module globals
