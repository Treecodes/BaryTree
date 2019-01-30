OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o \
	 treedriver.o partition.o xmalloc.o sort.o 
TARGET = tree.exe

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct.exe

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin.exe

#OFILESLIB = treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-batch.o partition.o xmalloc.o sort.o 
#TARGET = tree.exe

OUTDIR = ../bin/
CC = mpicc
#CC = /usr/local/Cellar/open-mpi/4.0.0/bin/mpicc

UNAME := $(shell uname -s)

ifeq ($(UNAME),Linux)
CFLAGS = -fopenmp -O3 -lm -std=gnu99 -fPIC
LDFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5
endif

ifeq ($(UNAME),Darwin)
MKLROOT = /opt/intel/mkl
CFLAGS = -Wall -O3 -fopenmp -lm -I${MKLROOT}/include -std=gnu99
LDFLAGS =  -DMKL_ILP64 -m64 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5
#LDFLAGS =  -DMKL_ILP64 -m64 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_core -ldl
#LDFLAGS =  -DMKL_ILP64 -m64 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_core -ldl -lpthread

#CFLAGS = -Wall -fopenmp -std=gnu99 -DMKL_ILP64 -m64 -I${MKLROOT}/include
#LDFLAGS = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif


#note: if mpicc doesn't properly pull the libraries and headers, look here:
#CFLAGS = -I/usr/local/Cellar/open-mpi/4.0.0/include
#LDFLAGS = -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/4.0.0/lib -lmpi

all: $(TARGET) $(TARGETA) $(TARGETB) bin library

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILESA) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(CFLAGS) $(OFILESB) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) $(OFILESB) $(TARGETB) 
	rm -f $(OUTDIR)*
	rm -f ../lib/*

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
