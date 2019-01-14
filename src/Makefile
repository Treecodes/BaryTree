OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-batch.o \
	 treedriver.o partition.o xmalloc.o sort.o 
TARGET = tree.exe

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct.exe

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin.exe

OUTDIR = ../bin/
CC = mpicc
#CC = /usr/local/Cellar/open-mpi/4.0.0/bin/mpicc

UNAME := $(shell uname -s)

ifeq ($(UNAME),Linux)
#CFLAGS = -acc -Minfo=accel -ta=tesla,time 
CFLAGS = -fopenmp
LIBS = -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

ifeq ($(UNAME),Darwin)
CFLAGS = -Wall -fopenmp
MKLROOT = /opt/intel/mkl
LIBS = -L/opt/intel/lib -L/opt/intel/mkl/lib -liomp5 -lmkl_intel_thread -lmkl_intel_ilp64 -lmkl_core -lpthread -lm -ldl
endif

#LIBS = ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl


all: $(TARGET) $(TARGETA) $(TARGETB) bin

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(OFILESA) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(CFLAGS) $(OFILESB) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) $(OFILESB) $(TARGETB) 
	rm -f $(OUTDIR)*

bin:                                                                                         
	cp $(TARGET) $(OUTDIR)$(TARGET)
	cp $(TARGETA) $(OUTDIR)$(TARGETA)
	cp $(TARGETB) $(OUTDIR)$(TARGETB)

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
tree-pc-yuk.o: tree-pc-yuk.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
xmalloc.o: xmalloc.c xmalloc.h
sort.o: sort.c
main-direct.o: main-direct.c array.h xmalloc.h tools.h
txt2bin.o: txt2bin.c array.h xmalloc.h
