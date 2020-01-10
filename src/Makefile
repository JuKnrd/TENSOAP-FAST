.PHONY: all clean

FF=gfortran
FOPTS=-O3 -g -fopenmp
LD=-lblas -llapack
OBJECTS1=sagpr_get_kernel.o
OBJECTS2=sagpr_train.o

all: mkbin sagpr_get_kernel sagpr_train

mkbin:
	if [ ! -d ../bin ];then mkdir ../bin;fi

sagpr_get_kernel: $(OBJECTS1)
	$(FF) $(FOPTS) -o ../bin/sagpr_get_kernel $(OBJECTS1) $(LD)

sagpr_train: $(OBJECTS2)
	$(FF) $(FOPTS) -o ../bin/sagpr_train $(OBJECTS2) $(LD)

%.o: %.f90 
	$(FF) $(FOPTS) -c $< 

clean: 
	rm -f *.o *.mod ../bin/*;rmdir ../bin
