FC=gfortran
#Note that command line reading does not seem to work with gcc 9. No idea why. 8 and 11 ok
FOPT= -O3 -fPIC #-std=legacy# -Wall -fbounds-check -g
# FOPT = -pg
SRCDIR = ./src
# NUMDIR = ./numerical
# QAGDIR = ./numerical/dqag
# TSDIR = ./numerical/TSPACK
# WDIR = ./Wfunctions
# RDIR = ./Rfunctions

MAIN = cosmion.o
MFSHR = init_conds.o star.o walk.o
NUMFOBJ =  num.o  rkf45.o rkf45full.o rkf45fullhistory.o newton.o #rkf45full_n2.o



# TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o


csharedlib.so: $(MFSHR)  $(NUMFOBJ)
	$(FC) $(FOPT) -shared -o $@ $(NUMFOBJ) $(MFSHR)

cosmion.x: $(MAIN)  csharedlib.so
	${FC} $(FOPT) -L. -Wl,-rpath,. -o cosmion.x -fopenmp $(MAIN) csharedlib.so


# debug: $(NUMFOBJ) $(MFSHR) $(MAIN)
# 		${FC} $(FOPT)  -Wl,-rpath,. -o cosmionDB.x $(MAIN) $(MFSHR) $(NUMFOBJ)



$(NUMFOBJ): %.o : $(SRCDIR)/%.f90
	$(FC) $(FOPT) -c  $<

$(MFSHR): %.o : $(SRCDIR)/%.f90
	$(FC) $(FOPT) -c  $< -fopenmp

$(MAIN): %.o : $(SRCDIR)/%.f90
	$(FC) $(FOPT) -c  $<

clean:
	rm -f *.o *.mod *.so cosmion.x
