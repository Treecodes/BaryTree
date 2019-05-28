#OFILES = main.o tools.o tree-cp.o tree-pc.o tree-batch.o \
#	 treedriver.o partition.o xmalloc.o sort.o
OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-hermite.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o \
	 treedriver.o partition.o xmalloc.o sort.o 
TARGET = tree.exe

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct.exe

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin.exe

#OFILESLIB = treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-batch.o partition.o xmalloc.o sort.o 
#TARGET = tree.exe

OUTDIR = ../bin/
#CC = /usr/local/Cellar/open-mpi/4.0.0/bin/mpicc
  
UNAME := $(shell uname -s)

#ifeq ($(UNAME),Linux)
#CC = mpicc
#CFLAGS = -fopenmp -O3 -g -lm -std=gnu99 -fPIC  
#LDFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5
#endif

ifeq ($(UNAME),Linux)
CC = mpicc -fPIC
#CFLAGS = -O3 -acc=strict -ta=nvidia,time -Minfo=accel -lm -mp  # to include openACC details (during compilation and runtime)
#CFLAGS = -O3 -acc -ta=nvidia,time -Minfo -lm -mp  # to include openACC details (during compilation and runtime)
#CFLAGS = -O3 -acc -ta=nvidia,time,keepptx,keepgpu -Minfo=accel -lm -mp  # to include openACC details (during compilation and runtime)

CFLAGS = -O3 -acc -ta=nvidia -lm -mp  # to supress any openACC details
CFLAGS35 = -O3 -acc -ta=tesla:cc35 -lm -mp  # to supress any openACC details
CFLAGS70 = -O3 -acc -ta=tesla:cc70 -lm -mp  # to supress any openACC details
#CFLAGS = -O3 -acc -ta=tesla:cc70 -lm -mp  # to compile for the titanv

#CFLAGS = -O3 -ta=nvidia,time -Minfo=accel -lm  # to supress any openACC details and not use pragmas
CXXFLAGS='-fPIC'

 
## While testing openACC, drop the mkl business
#CFLAGS = -fopenmp -O3 -lm -I${MKLROOT}/include -std=gnu99
#LDFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5
endif

ifeq ($(UNAME),Darwin)
CC = mpicc
MKLROOT = /opt/intel/mkl
CFLAGS = -Wall -O3 -fopenmp -lm -I${MKLROOT}/include -std=gnu99
LDFLAGS =  -DMKL_ILP64 -m64 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5
endif 
 

#note: if mpicc doesn't properly pull the libraries and headers, look here:
#CFLAGS = -I/usr/local/Cellar/open-mpi/4.0.0/include
#LDFLAGS = -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/4.0.0/lib -lmpi

#cc35Lib: $(TARGET) $(TARGETA) $(TARGETB) library35

	
#library35: treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o partition.o xmalloc.o sort.o
#	$(CC35) $(CFLAGS) $(LDFLAGS) treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o partition.o xmalloc.o sort.o -shared -o ../lib35/libtreedriverWrapper.so


#all: $(TARGET) $(TARGETA) $(TARGETB) library bin
all: $(TARGET) $(TARGETA) $(TARGETB) bin

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILESA) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(CFLAGS) $(OFILESB) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) $(OFILESB) $(TARGETB) #$(library)
	rm -f $(OUTDIR)*
	rm -f ../libTest/*
	rm -f *.so    
	rm -f *.o 

bin:                                                                                            
	mv $(TARGET) $(OUTDIR)$(TARGET)
	mv $(TARGETA) $(OUTDIR)$(TARGETA)
	mv $(TARGETB) $(OUTDIR)$(TARGETB)
	

library: treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o partition.o xmalloc.o sort.o
	$(CC) $(CFLAGS) $(LDFLAGS) treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o partition.o xmalloc.o sort.o -shared -o ../lib/libtreedriverWrapper.so

main.o: main.c array.h xmalloc.h treedriver.h tools.h particles.h
tools.o: tools.c array.h xmalloc.h
treedriver.o: treedriver.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h tree.h treedriver.h
partition.o: partition.c tools.h partition.h
tree-cp.o: tree-cp.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-batch.o: tree-batch.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h array.h
tree-cp-yuk.o: tree-cp-yuk.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc.o: tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-hermite.o: tree-pc-hermite.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-yuk.o: tree-pc-yuk.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-yuk-SS.o: tree-pc-yuk-SS.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-coulomb-SS.o: tree-pc-coulomb-SS.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
xmalloc.o: xmalloc.c xmalloc.h
sort.o: sort.c
main-direct.o: main-direct.c array.h xmalloc.h tools.h
txt2bin.o: txt2bin.c array.h xmalloc.h
treedriverWrapper.o: treedriverWrapper.c
