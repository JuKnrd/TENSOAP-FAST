.PHONY: all clean

FF=gfortran
FOPTS=-O3 -g -fopenmp
OBJECTS1=sagpr_get_kernel.o

all: mkbin sagpr_get_kernel

mkbin:
	if [ ! -d ../bin ];then mkdir ../bin;fi

sagpr_get_kernel: $(OBJECTS1)
	$(FF) $(FOPTS) -o ../bin/sagpr_get_kernel $(OBJECTS1)

%.o: %.f90 
	$(FF) $(FOPTS) -c $< 

clean: 
	rm -f *.o *.mod ../bin/*;rmdir ../bin
