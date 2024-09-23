module globals
   use, intrinsic :: iso_c_binding
   implicit none
!=======================================================================
!  ARRAY DIMENSIONS:
   integer, parameter :: maxeq=100000           !maximum number of equations
   integer, parameter :: maxz=30                !maximum number of axial nodes
   integer, parameter :: maxr=30                !maximum number of radial nodes
   integer, parameter :: maxtab=50              !maximum number of entries in tables

!=======================================================================
!  DIMENSIONS OF TABLES:
   integer nqv                                  !number of entries in table of fuel rod power density vs time
   integer nfsto                                !number of fuel stoichiometry vs fuel burnup table entries
   integer ntco                                 !number of fuel rod clad outer temperature table entries
   integer ntple                                !number of gas plenum temperature table entries

!=======================================================================
!  FLAGS:
   integer ifgr                                 !option of fission gas release calculation (0-no, 1-yes)
   integer ifcreep                              !option of fuel creep calculation (0-no, 1-yes)
   integer ifswel                               !option of fuel swelling calculation (0-no, 1-yes)
   integer iccreep                              !option of clad creep calculation (0-no, 1-yes)
   integer ifreloc                              !option of fuel relocation calculation (0-no, 1-yes)
   integer inomech                              !option of no mechanics calculation (0-no, 1-yes)

   character*4 flag(maxz)                       !flag of gap state ('clos', 'open')

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
   real*8 tend                                  !terminal time (s)
   real*8 dtout                                 !output time step (s)
   real*8 hmax                                  !maximum time step (s)

!=======================================================================
!  GEOMETRY PARAMETERS:
   integer nf                                   !number of radial fuel nodes
   integer nc                                   !number of radial clad nodes
   integer nzz                                  !number of axial layers

   real*8 gap(maxz)                             !radial fuel-clad gap (m)
   real*8 gapth(maxz)                           !"thermal" gap width (m)
   real*8 reloc(maxz)                           !fraction of the gap closed due to fuel relocation
   real*8 gask(maxz)                            !gas thermal conductivity (W/mK)
   real*8 rad0(maxr)                            !initial node radius (m)
   real*8 rad_0(maxr)                           !initial radius of boundary between nodes(m)
   real*8 rad(maxr,maxz)                        !node radius (m)
   real*8 rfi0                                  !initial inner fuel radius (m)
   real*8 rfi(maxz)                             !inner fuel radius (m)
   real*8 rfo0                                  !initial outer fuel radius (m)
   real*8 rfo(maxz)                             !outer fuel radius (m)
   real*8 rci0                                  !initial inner clad radius (m)
   real*8 rci(maxz)                             !inner clad radius (m)
   real*8 rco0                                  !initial outer clad radius (m)
   real*8 rco(maxz)                             !outer clad radius (m)
   real*8 dz0(maxz)                             !initial node height (m)
   real*8 dzf(maxz)                             !fuel node height (m)
   real*8 dzc(maxz)                             !clad node height (m)
   real*8 drf0                                  !initial fuel node thickness (m)
   real*8 drf(maxr,maxz)                        !fuel node thickness (m)
   real*8 drc0                                  !initial clad node thickness (m)
   real*8 drc(maxr,maxz)                        !clad node thickness (m)
   real*8 az0(maxr)                             !initial node cross sectional area (m2)
   real*8 azf0                                  !initial fuel cross sectional area (m2)
   real*8 azc0                                  !initial clad cross sectional area (m2)
   real*8 vgp                                   !gas plenum volume (m3)
   real*8 h0                                    !initial fuel height (m)

   character*6 fmat
   character*6 gmat
   character*6 cmat

!=======================================================================
!  TEMPERATURE PARAMETERS:
   real*8 tem0                                  !initial temperature (K)
   real*8 tem(maxr,maxz)                        !temperature in node (K)
   real*8 dtem(maxr,maxz)                       !time derivative of temperature in node (K)
   real*8 tco(maxz,maxtab)                      !fuel rod clad outer temperature (K) in table of fuel rod clad outer temperature vs time
   real*8 ttco(maxtab)                          !time (s) in table of fuel rod clad outer temperature vs time

!=======================================================================
!  ENERGY PARAMETERS:
   real*8 qv(maxz,maxtab)                       !fuel rod power density (W/m3) in table of fuel rod power density vs time
   real*8 tqv(maxtab)                           !time (s) in table of fuel rod power density vs time
   real*8 ql(maxz)                              !fuel rod linear power (W/m)
   real*8 qqv1(maxz)                            !fuel power density (W/m3)

!=======================================================================
!  GAS PLENUM PARAMETERS:
   real*8 gpres0                                !initial inner gas pressure (Pa)
   real*8 gpres                                 !inner gas pressure (Pa)
   real*8 mu0                                   !as fabricated gas amount (mol)
   real*8 tple(maxtab)                          !gas plenum temperature (K) in table of gas plenum temperature vs time
   real*8 ttple(maxtab)                         !time (s) in table of gas plenum temperature vs time

!=======================================================================
!  COOLANT PARAMETERS (FOR BASE IRRADIATION ANALYSIS):
   real*8 pcool                                 !coolant pressure (Pa)

!=======================================================================
!  MATERIAL PROPERTIES PARAMETERS:
   real*8 bup(maxz)                             !fuel burnup (MWd/kgU)
   real*8 dbup(maxz)                            !time derivative of fuel burnup (MWd/kgU/s)
   real*8 pucont                                !plutonium content in fuel (-)
   real*8 sto0                                  !initial fuel stoichiometry (-)
   real*8 rof0                                  !initial fuel density (kg/m3)
   real*8 roc0                                  !initial clad density (kg/m3)
   real*8 rof(maxr,maxz)                        !fuel density (kg/m3)
   real*8 kfuel(maxr,maxz)                      !fuel thermal condcutivity (W/mK)
   real*8 kclad(maxr,maxz)                      !clad thermal condcutivity (W/mK)
   real*8 bsto(maxtab)                          !burnup (MWd/kgHM) in table of fuel stoichiometry vs. burnup
   real*8 stob(maxtab)                          !fuel stoichiometry (-) in table of fuel stoichiometry vs. burnup

!=======================================================================
!  GAP PARAMETERS:
   real*8 pfc(maxz)                             !pellet-to-cladding contact pressure (MPa)
   real*8 dpfc(maxz)                            !time derivative of pellet-to-cladding contact pressure (MPa/s)
   real*8 rufc                                  !cladding inner surface roughness (m)
   real*8 ruff                                  !fuel outer surface roughness (m)
   real*8 hgapt(maxz)                           !gap conductance (W/m**2-K)
   real*8 hgapi(3,maxz)                         !components of gap conductance (W/m**2-K)
   real*8 xmol(5)                               !number of moles of gas components (-)
   real*8 ajump(maxz)                           !solid-gas temperature jump (m)

!=======================================================================
!  CONSTANTS:
   real*8, parameter :: rmu=8.314d0             !universal gas constant
   real*8, parameter :: avo=6.0247d23           !Avogadro number
   real*8, parameter :: pi=3.141592653589793d0  !pi number

!=======================================================================
!  MECHANICAL PARAMETERS:
   real*8 sigh(maxr,maxz)                       !clad hoop stress (Pa)
   real*8 sigr(maxr,maxz)                       !clad radial stress (Pa)
   real*8 sigz(maxr,maxz)                       !clad axial stress (Pa)
   real*8 sig(maxr,maxz)                        !clad effective stress (Pa)
   real*8 eh(maxr,maxz)                         !clad total hoop strain (m/m)
   real*8 er(maxr,maxz)                         !clad total radial strain (m/m)
   real*8 ez(maxz)                              !clad total axial strain (m/m)
   real*8 deh(maxr,maxz)                        !clad total hoop strain time derivative (m/m)
   real*8 dez(maxz)                             !clad total axial strain time derivative (1/s)
   real*8 ece(maxr,maxz)                        !clad creep effective strain (m/m)
   real*8 ech(maxr,maxz)                        !clad creep hoop strain (m/m)
   real*8 ecr(maxr,maxz)                        !clad creep radial strain (m/m)
   real*8 ecz(maxr,maxz)                        !clad creep axial strain (m/m)
   real*8 dece(maxr,maxz)                       !time derivative of clad creep effective strain (m/m)
   real*8 dech(maxr,maxz)                       !time derivative of clad creep hoop strain (m/m)
   real*8 decr(maxr,maxz)                       !time derivative of clad creep radial strain (m/m)
   real*8 decz(maxr,maxz)                       !time derivative of clad creep axial strain (m/m)
   real*8 et(maxr,maxz)                         !clad thermal linear strain (m/m)
   real*8 sigfh(maxr,maxz)                      !fuel hoop stress (Pa)
   real*8 sigfr(maxr,maxz)                      !fuel radial stress (Pa)
   real*8 sigfz(maxr,maxz)                      !fuel axial stress (Pa)
   real*8 sigf(maxr,maxz)                       !fuel effective stress (Pa)
   real*8 efh(maxr,maxz)                        !fuel total hoop strain (m/m)
   real*8 efr(maxr,maxz)                        !fuel total radial strain (m/m)
   real*8 efz(maxz)                             !fuel total axial strain (m/m)
   real*8 defh(maxr,maxz)                       !fuel total hoop strain time derivative (m/m)
   real*8 defz(maxz)                            !fuel total axial strain time derivative (1/s)
   real*8 eft(maxr,maxz)                        !fuel linear thermal expansion (m/m)
   real*8 efs(maxr,maxz)                        !fuel volumetric swelling strain (m/m)
   real*8 defs(maxr,maxz)                       !time derivative of fuel volumetric swelling strain (m/m)
   real*8 fcreepc(3)                            !user input constants for the fuel creep rate law
   real*8 efce(maxr,maxz)                       !fuel creep effective strain (m/m)
   real*8 efch(maxr,maxz)                       !fuel creep hoop strain (m/m)
   real*8 efcr(maxr,maxz)                       !fuel creep radial strain (m/m)
   real*8 efcz(maxr,maxz)                       !fuel creep axial strain (m/m)
   real*8 defce(maxr,maxz)                      !time derivative of fuel creep effective strain (1/s)
   real*8 defch(maxr,maxz)                      !time derivative of fuel creep hoop strain (1/s)
   real*8 defcr(maxr,maxz)                      !time derivative of fuel creep radial strain (1/s)
   real*8 defcz(maxr,maxz)                      !time derivative of fuel creep axial strain (1/s)

!=======================================================================
! FISSION GAS RELEASE AND SWELLING PARAMETERS:
   real*8 fggen                                 !amount of fission gas produced in fuel rod (mol/s)
   real*8 fgrel                                 !amount of fission gas released in fuel rod (mol)
   real*8 dfggen                                !time derivative of amount of fission gas produced in fuel rod (mol/s)
   real*8 dfgrel                                !time derivative of amount of fission gas released in fuel rod (mol/s)

end module globals
