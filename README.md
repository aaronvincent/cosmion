# cosmion

### Quickstart

If you have gfortran (installed with gcc), cosmion should work (almost) out of the box. In the main directory, type

`bash runandtest.sh`

This should compile the main program (cosmion) and run a test example. Have a look in the file to see the commands that it uses. To rerun the code once you have compiled, you can write e.g.

`./cosmion.x 5. 1.d-37 100000 "positions.dat" nucleonSD const`

In order the arguments are 

1) dark matter mass (GeV)
2) DM-nucleon cross section (cm^2)
3) Number of collisions
4) Output file name
5) nucleonSI (spin-independent) or nucleonSD (spin-dependent). electron is an option, but has not been fully tested yet. Use at your own risk
6) const (constant cross section), or v2, v4, q2, q4, if you want velocity or momentum-dependent interactions

The output will be stored in the positions.dat file. Columns in this file are: 

1) The radial coordinate of the particle at the moment of collision (cm). Simulations in an analytic, approximated gravitational potential use this column as the x-coordinate position of the particle (cm). Note that the centre of the star is set to (0,0,0) in Cartesian coordinates.
2) Simulations in an analytic, approximated gravitational potential use this column as the y-coordinate position of the particle (cm). This column is set to zero in a realistic potential.
3) Simulations in an analytic, approximated gravitational potential use this column as the z-coordinate position of the particle (cm). This column is set to zero in a realistic potential.
4) The radial velocity of the particle going into the collision (cm/s). Simulations in an analytic, approximated gravitational potential use this column as the x-coordinate velocity of the particle going into the collision (cm/s).
5) The angular velocity of the particle going into the collision (cm/s). Simulations in an analytic, approximated gravitational potential use this column as the y-coordinate velocity of the particle going into the collision (cm/s).
6) Simulations in an analytic, approximated gravitational potential use this column as the z-coordinate velocity of the particle going into the collision (cm/s). This column is set to zero in a realistic potential.
7) The x-coordinate velocity of the particle coming out of the collision (cm/s).
8) The y-coordinate velocity of the particle coming out of the collision (cm/s).
9) The z-coordinate velocity of the particle coming out of the collision (cm/s).
10) The time between this collision and the previous collision (s).
11) A flag indicating whether the particle has left the star. If the flag is 0, the particle has collided within the star. If the flag is 1, this row records the properties of the particle as enters the surface of the star. If the flag is 2, the particle has escaped the star's gravitational potential, and this row records the initial properties of a newly spawned particle. If the flag is 3, this row records the properties of the particle as exits the surface of the star.
12) The rate a which the particle will collide at this radius and speed (s^-1).

Note that the first row indicates the initial properties of the spawned particle.

The accompanying example_processing_script.py can be used to calculate the luminosity profile from particle collision data in the positions.dat output files. For further details, see the script's docstring.

Note: If you get GOMP-related errors, it's because you have not installed openmp. You can either do so, or comment out the `-fopenmp` option in the Makefile:

`	$(FC) $(FOPT) -c  $< # -fopenmp`
