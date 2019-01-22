OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o tree-batch.o \
	 treedriver.o partition.o xmalloc.o sort.o 
TARGET = tree.exe

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct.exe

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin.exe

OUTDIR = ../bin/
CC = mpicc 

UNAME := $(shell uname -s)

ifeq ($(UNAME),Linux)
CFLAGS = -acc -Minfo=accel -ta=tesla,time
endif
ifeq ($(UNAME),Darwin)
CFLAGS = -Wall
endif

#note: if mpicc doesn't properly pull the libraries and headers, look here:
#CFLAGS = -I/usr/local/Cellar/open-mpi/4.0.0/include
#LDFLAGS = -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/4.0.0/lib -lmpi

all: $(TARGET) $(TARGETA) $(TARGETB) bin

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILESA) -o $@

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
