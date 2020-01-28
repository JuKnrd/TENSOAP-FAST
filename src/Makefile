.PHONY: all clean

FF=gfortran
FOPTS=-O3 -g -fopenmp -Wunused -funroll-loops
LD=-lblas -llapack
OBJECTS1=mod_sagpr.o sagpr_get_PS.o
OBJECTS2=mod_sagpr.o sagpr_get_kernel.o
OBJECTS3=mod_sagpr.o sagpr_predict.o
OBJECTS4=mod_sagpr.o
OBJECTS5=sagpr_apply.o

all: mkbin sagpr_get_PS sagpr_get_kernel sagpr_predict libsoapfast sagpr_apply

mkbin:
	@if [ ! -d ../bin ];then mkdir ../bin;fi

sagpr_get_PS: $(OBJECTS1)
	$(FF) $(FOPTS) -o ../bin/sagpr_get_PS $(OBJECTS1) $(LD)

sagpr_get_kernel: $(OBJECTS2)
	$(FF) $(FOPTS) -o ../bin/sagpr_get_kernel $(OBJECTS2) $(LD)

sagpr_predict: $(OBJECTS3)
	$(FF) $(FOPTS) -o ../bin/sagpr_predict $(OBJECTS3) $(LD)

libsoapfast:
	@if [ ! -d ../lib ];then mkdir ../lib;fi
	ar r ../lib/libsoapfast.a $(OBJECTS4)

sagpr_apply: $(OBJECTS5)
	$(FF) $(FOPTS) -o ../bin/sagpr_apply $(OBJECTS5) -L../lib -lsoapfast $(LD)

%.o: %.f90 
	$(FF) $(FOPTS) -c $< 

install: all mkinstall

mkinstall:
	@if [ "${SAGPR_PREFIX}" != "" ];then cp ../bin/sagpr_* ${SAGPR_PREFIX}/;fi

clean: 
	rm -rf *.o *.mod ../bin/* ../bin ../lib
