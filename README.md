# cosmion

### Quickstart

If you have a working modern Fortran compiler, cosmion should work (almost) out of the box. In the main directory, type

`bash runandcrash.sh`

This should compile the main program (cosmion) and run a test example. Have a look in the file to see the commands that it uses. To rerun the code once you have compiled, you can write e.g.

`./cosmion.x 5. 1.d-37 100000 "positions.dat" SD'

In order the arguments are 

1) dark matter mass (GeV)
2) DM-nucleon cross section (cm^2)
3) Number of collisions
4) Output file name
5) SI (Spin-independent) or SD (spin-dependent)

The output will be stored in the positions.dat file. Columns in this file are: 


Note: If you get GOMP-related errors, it's because you have not installed openmp. You can either do so, or comment out the `-fopenmp` option in the Makefile:

`	$(FC) $(FOPT) -c  $< # -fopenmp`
