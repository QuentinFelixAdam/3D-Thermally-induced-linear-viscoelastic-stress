The presented code allows the calculations of three-dimensional thermally induced thermo-viscoelastic stress. It was utilized to quantify thermal stress at the surface of an asphalt pavement. Asphalt was viewed as a linear viscoelastic (solid), and a constitutive relation was selected for this purpose, governed by a constant Poisson's ratio ("nu" in the code, unitless), and a thermo-sensitive relaxation modulus associated with a certain temperature level. Here, Delta_T = T - Tr (degrees Celsius) wherein T is the surface temperature and Tr serves as a reference temperature for the thermal analysis. The meaning of Tr is as follows: if the material is exposed to a constant temperature level above Tr it can approach a stress-free state provided sufficient relaxation time is allowed. Conversely, if the material is exposed to a constant temperature level below Tr, induced stresses do not approach zero regardless of the amount of relaxation time allowed.

The following strategy is adopted: stress calculations are commenced from a stress-free state whenever T < Tr, i.e. the surface temperature drops below the reference temperature. As long as the surface temperature remained below Tr, calculation continues. If, at a certain time instance T > Tr, i.e. the surface temperature climbs above the reference temperature, and also thermal stress is tensile, the calculation continued. However, if at a certain time instance T > Tr (i.e. the surface temperature climbs above the reference temperature) and also thermal stress is compressive, the stress calculation was stopped -- only to be resumed later from a stress-free state when the surface temperature again drops below the reference temperature. The time-temperature superposition principle was applied by assuming thermo-rheological simplicity of the asphalt. The relaxation function associated with the linear viscoelastic solid taken as a four-parameter sigmoid (in double log scale). Einf (MPa) is the long-term equilibrium modulus, TauD (s) and nD (unitless) are shape parameters, and E0 (MPa) is the instantaneous modulus. All required parameters are embedded within the code.

The code, named "ThermalVE.f90", is implemented in Fortran and takes as input a file "Tsurface.dat" which contains values arranged in column (Unix LF format). Those values represent temperature levels. It is assumed that at least two values are contained in "Tsurface.dat". The output is contained within the file "VEStress.dat", in which the thermal stress values are stored. Every execution of "ThermalVE.f90" writes one value in "VEStress.dat", which corresponds to the thermal stress induced by the last temperature reading in "Tsurface.dat".

The compiler is GNU Fortran (GCC) Version 4.8.5 20150623. There is a possibility to parallelize the code with OpenMP.

In order to run the code, make sure all files are in the same folder ("ThermalVE.f90" and "Tsurface.dat"), open a terminal to navigate to the folder in which the files are placed, and type "f90 -free ThermalVE.f90 -o exe" followed by "exe". An example file "Tsurface.dat" is provided.

The code is validated against measurements. Please see the file "Thermally-induced viscoelastic stresses calculation validation.txt".

This code is based on a research paper to be published any time soon -- hopefully. Should you want to reuse that code and this paper be published, I will embed here a citation.

For any question: contact me via quentin.f.adam@gmail.com
