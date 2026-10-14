#OPT += -DEXACT_HERNQUIST   #generates a Hernquist sphere based on distribution function (no approximation)
#OPT += -DOUTPUT_DF         #writes out dM/dE and f(E) and the exact values for comparison (see DF.txt)

OPTIONS = $(OPT)

EXEC   = MakeDiskGalaxy

OBJS   = main.o cooling.o read_param.o effmodel.o  globvars.o \
	halo.o haloset.o escape.o \
	disk.o diskset.o \
	bulge.o bulgeset.o \
	gas.o  gasset.o \
	save.o init.o structure.o \
	misc.o toomre.o force.o \
	forcetree.o \
	nrsrc/bessi0.o   nrsrc/bessj1.o   \
	nrsrc/erff.o     nrsrc/polint.o   nrsrc/spline.o \
	nrsrc/bessi1.o   nrsrc/bessk0.o   nrsrc/gasdev.o  \
	nrsrc/qromb.o    nrsrc/splint.o  \
	nrsrc/bessj0.o   nrsrc/bessk1.o   \
	nrsrc/nrutil.o   nrsrc/ran1.o     nrsrc/trapzd.o \
	nrsrc/gammp.o    nrsrc/gcf.o      nrsrc/gser.o  \
        nrsrc/gammln.o   nrsrc/indexx.o   nrsrc/zriddr.o


INCL   = prototypes.h globvars.h nrsrc/nrutil.h nrsrc/nr.h


CFLAGS =  $(OPTIONS) -O2 -Wall -g

LIBS   =  -lm  -g

CC     =  cc


$(EXEC): $(OBJS) 
	cc $(CFLAGS) $(OBJS) -o $(EXEC)  $(LIBS)
	mv MakeDiskGalaxy MakeDiskGalaxy_ext

$(OBJS): $(INCL) 

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
	rm MakeDiskGalaxy_ext











