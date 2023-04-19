# cosmion

### Quickstart

If you have a working modern Fortran compiler, cosmion should work (almost) out of the box. In the main directory, type

`.bash runandcrash.sh`

This should compile the main program (cosmion) and run a test example. Have a look in the file to see the commands that it uses. If you get GOMP-related errors, it's because you have not installed openmp. You can either do so, or comment out the `-fopenmp` line in the Makefile:

`$(MFSHR): %.o : $(SRCDIR)/%.f90`
`	$(FC) $(FOPT) -c  $< # -fopenmp`
