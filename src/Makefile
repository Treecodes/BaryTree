#OFILES = main.o tools.o tree-cp.o tree-pc.o tree-batch.o \
#	 treedriver.o partition.o xmalloc.o sort.o
OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-hermite.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o \
	 treedriver.o partition.o xmalloc.o sort.o tree-pc-hermite-coulombSS.o
TARGET = tree.exe

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct.exe

OFILESC = main-direct-distributed.o tools.o xmalloc.o
TARGETC = direct-distributed.exe

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin.exe

OUTDIR = ../bin/


OFILES_cpu = main_cpu.o tools_cpu.o tree-cp_cpu.o tree-cp-yuk_cpu.o tree-pc_cpu.o tree-pc-hermite_cpu.o tree-pc-yuk_cpu.o tree-pc-yuk-SS_cpu.o \
	tree-pc-coulomb-SS_cpu.o tree-batch_cpu.o \
	treedriver_cpu.o partition_cpu.o xmalloc_cpu.o sort_cpu.o tree-pc-hermite-coulombSS_cpu.o
TARGET_cpu = tree_cpu.exe

OFILESA_cpu = main-direct_cpu.o tools_cpu.o xmalloc_cpu.o
TARGETA_cpu = direct_cpu.exe

OFILESC_cpu = main-direct-distributed_cpu.o tools_cpu.o xmalloc_cpu.o
TARGETC_cpu = direct-distributed_cpu.exe

OFILESB_cpu = txt2bin_cpu.o sort_cpu.o xmalloc_cpu.o
TARGETB_cpu = txt2bin_cpu.exe

OUTDIR_cpu = ../bin_noACC/
#CC = /usr/local/Cellar/open-mpi/4.0.0/bin/mpicc
  
UNAME := $(shell uname -s)



ifeq ($(UNAME),Linux)
CC_cpu = mpicc
CFLAGS_cpu = -fopenmp -g -lm -std=gnu99 -fPIC  
LDFLAGS_cpu = -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -ldl -liomp5

CC = mpicc -fPIC
#CFLAGS = -O3 -acc=strict -ta=nvidia,time -Minfo=accel -lm -mp  # to include openACC details (during compilation and runtime)
#CFLAGS = -O3 -acc -ta=nvidia,time -Minfo -lm -mp  # to include openACC details (during compilation and runtime)
#CFLAGS = -O3 -acc -ta=nvidia,time,keepptx,keepgpu -Minfo=accel -lm -mp  # to include openACC details (during compilation and runtime)

CFLAGS = -O3 -acc -ta=nvidia -lm -mp  # to supress any openACC details
#CFLAGS35 = -O3 -acc -ta=tesla:cc35 -lm -mp  # to supress any openACC details
#CFLAGS70 = -O3 -acc -ta=tesla:cc70 -lm -mp  # to supress any openACC details
#CFLAGS = -O3 -acc -ta=tesla:cc70 -lm -mp  # to compile for the titanv

#CFLAGS = -O3 -ta=nvidia,time -Minfo=accel -lm  # to supress any openACC details and not use pragmas
CXXFLAGS='-fPIC'

# While testing openACC, drop the mkl business
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
gpu: $(TARGET) $(TARGETA) $(TARGETB) $(TARGETC) bin library

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILESA) -o $@
	
$(TARGETC): $(OFILESC)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILESC) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(CFLAGS) $(OFILESB) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) $(OFILESB) $(TARGETB) #$(library)
	rm -f $(OFILES_cpu) $(TARGET_cpu) $(OFILESA_cpu) $(TARGETA_cpu) $(OFILESB_cpu) $(TARGETB_cpu) #$(library)
	rm -f $(OUTDIR)*
	rm -f $(OUTDIR_cpu)*
	rm -f ../libTest/*
	rm -f *.so    
	rm -f *.o 

bin:                                                                                            
	mv $(TARGET) $(OUTDIR)$(TARGET)
	mv $(TARGETA) $(OUTDIR)$(TARGETA)
	mv $(TARGETB) $(OUTDIR)$(TARGETB)
	mv $(TARGETC) $(OUTDIR)$(TARGETC)
	
cpu: $(TARGET_cpu) $(TARGETA_cpu) $(TARGETB_cpu) $(TARGETC_cpu) bin_cpu library_cpu

$(TARGET_cpu): $(OFILES) 
	$(CC_cpu) $(CFLAGS_cpu) $(LDFLAGS_cpu) $(OFILES_cpu) -o $@

$(TARGETA_cpu): $(OFILESA_cpu)
	$(CC_cpu) $(CFLAGS_cpu) $(LDFLAGS_cpu) $(OFILESA_cpu) -o $@
	
$(TARGETC_cpu): $(OFILESC_cpu)
	$(CC_cpu) $(CFLAGS_cpu) $(LDFLAGS_cpu) $(OFILESC_cpu) -o $@

$(TARGETB_cpu): $(OFILESB_cpu)
	$(CC_cpu) $(CFLAGS_cpu) $(OFILESB_cpu) -o $@
	
bin_cpu:                                                                                            
	mv $(TARGET_cpu) $(OUTDIR_cpu)$(TARGET_cpu)
	mv $(TARGETA_cpu) $(OUTDIR_cpu)$(TARGETA_cpu)
	mv $(TARGETB_cpu) $(OUTDIR_cpu)$(TARGETB_cpu)
	mv $(TARGETC_cpu) $(OUTDIR_cpu)$(TARGETC_cpu)
	

library: treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-batch.o partition.o xmalloc.o sort.o
	$(CC) $(CFLAGS) $(LDFLAGS) treedriverWrapper.o treedriver.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-pc-yuk-SS.o tree-pc-coulomb-SS.o tree-pc-hermite.o tree-pc-hermite-coulombSS.o tree-batch.o partition.o xmalloc.o sort.o -shared -o ../lib/libtreedriverWrapper.so

library_cpu: treedriverWrapper_cpu.o treedriver_cpu.o tools_cpu.o tree-cp_cpu.o tree-cp-yuk_cpu.o tree-pc_cpu_cpu.o tree-pc-yuk_cpu.o tree-pc-yuk-SS_cpu.o tree-pc-coulomb-SS_cpu.o tree-batch_cpu.o partition_cpu.o xmalloc_cpu.o sort_cpu.o
	$(CC_cpu) $(CFLAGS_cpu) $(LDFLAGS_cpu) treedriverWrapper_cpu.o treedriver_cpu.o tools_cpu.o tree-cp_cpu.o tree-cp-yuk_cpu.o tree-pc_cpu.o tree-pc-yuk_cpu.o tree-pc-yuk-SS_cpu.o tree-pc-coulomb-SS_cpu.o tree-pc-hermite_cpu.o tree-pc-hermite-coulombSS_cpu.o tree-batch_cpu.o partition_cpu.o xmalloc_cpu.o sort_cpu.o -shared -o ../lib_cpu/libtreedriverWrapper.so


main_cpu.o: main.c array.h xmalloc.h treedriver.h tools.h particles.h
tools_cpu.o: tools.c array.h xmalloc.h
treedriver_cpu.o: treedriver.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h tree.h treedriver.h
partition_cpu.o: partition.c tools.h partition.h
tree-cp_cpu.o: tree-cp.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-batch_cpu.o: tree-batch.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h array.h
tree-cp-yuk_cpu.o: tree-cp-yuk.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc_cpu.o: tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-hermite_cpu.o: tree-pc-hermite.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-hermite-coulombSS_cpu.o: tree-pc-hermite-coulombSS.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-yuk_cpu.o: tree-pc-yuk.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-yuk-SS_cpu.o: tree-pc-yuk-SS.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
tree-pc-coulomb-SS_cpu.o: tree-pc-coulomb-SS.c tree-pc.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
		  tools.h partition.h tree.h
xmalloc_cpu.o: xmalloc.c xmalloc.h
sort_cpu.o: sort.c
main-direct_cpu.o: main-direct.c array.h xmalloc.h tools.h
main-direct-distributed_cpu.o: main-direct-distributed.c array.h xmalloc.h tools.h
txt2bin_cpu.o: txt2bin.c array.h xmalloc.h
treedriverWrapper_cpu.o: treedriverWrapper.c


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
tree-pc-hermite-coulombSS.o: tree-pc-hermite-coulombSS.c array.h xmalloc.h globvars.h tnode.h batch.h particles.h \
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
main-direct-distributed.o: main-direct-distributed.c array.h xmalloc.h tools.h
txt2bin.o: txt2bin.c array.h xmalloc.h
treedriverWrapper.o: treedriverWrapper.c



