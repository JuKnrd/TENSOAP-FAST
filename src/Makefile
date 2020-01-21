.PHONY: all clean

FF=gfortran
FOPTS=-O3 -g -fopenmp -Wunused
LD=#-lblas -llapack
OBJECTS1=mod_sagpr.o sagpr_get_PS.o
OBJECTS2=mod_sagpr.o sagpr_get_kernel.o
OBJECTS3=mod_sagpr.o sagpr_predict.o
OBJECTS4=mod_sagpr.o

all: mkbin sagpr_get_PS sagpr_get_kernel sagpr_predict libsagpr sagpr_apply

mkbin:
	@if [ ! -d ../bin ];then mkdir ../bin;fi

sagpr_get_PS:
	@echo "NOT YET IMPLEMENTED"

sagpr_get_kernel: $(OBJECTS2)
	$(FF) $(FOPTS) -o ../bin/sagpr_get_kernel $(OBJECTS2) $(LD)

sagpr_predict: $(OBJECTS3)
	$(FF) $(FOPTS) -o ../bin/sagpr_predict $(OBJECTS3) $(LD)

libsagpr:
	@if [ ! -d ../lib ];then mkdir ../lib;fi
	ar r ../lib/libsoapfast.a $(OBJECTS4)

sagpr_apply:
	@echo "NOT YET IMPLEMENTED"

%.o: %.f90 
	$(FF) $(FOPTS) -c $< 

clean: 
	rm -rf *.o *.mod ../bin/* ../bin ../lib
