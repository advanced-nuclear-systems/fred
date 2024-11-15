# FRED
## Fast reactor fuel behaviour code

Open-source Fortran code uses the SUNDIALS library for numerical modelling of base irradiation of a fast reactor fuel pin, accounting for fuel and clad heat transfer, fuel and clad stress-strain conditions, fuel and clad thermal expansion and creep, fuel swelling and fission gas release, evolution of inner gas pressure and composition as well as evolution of fuel-clad gap conductance. 
Algebraic and ordinary differential equations are solved by finite-difference method on a structured r-z cylindrical mesh using the SUNDIALS library.

As the input, FRED takes time-dependent axial profile of power density assuming flat power distribution over the radius, time-dependent axial profile of clad outer temperature, as-manufactured fuel and clad dimensions, mesh specification, as well as several constants for material properties and flags to activate specific models.

The typical dataset of the calculational results includes time-dependent radial and axial maps of fuel and clad deformations, stresses, and temperatures; axial profiles of fuel-clad gap conductance, contact pressure, burnup, fuel and clad dimensions.

Available materials: 
- Fuel: MOX.
- Clad: AIM1 stainless steel, T91 stainless steel.
- Bonding: helium, xenon, krypton.
- Inner gas: helium, xenon, krypton.

