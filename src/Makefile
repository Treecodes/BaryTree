OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o \
	 treedriver.o partition.o xmalloc.o sort.o
TARGET = tree

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin

OFILESC = makegrid.o
TARGETC = makegrid

CC = mpicc
CFLAGS = -Wall -O2 -std=c99
$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(OFILESA) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(OFILESB) -o $@

$(TARGETC): $(OFILESC)
	$(CC) $(OFILESC) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA)
	rm -f $(OFILESB) $(TARGETB) $(TARGETC) $(OFILESC)


main.o: main.c array.h xmalloc.h treedriver.h tools.h
tools.o: tools.c array.h xmalloc.h
treedriver.o: treedriver.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  tree.h treedriver.h
partition.o: partition.c tools.h partition.h
tree-cp.o: tree-cp.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
tree-cp-yuk.o: tree-cp-yuk.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
tree-pc.o: tree-pc.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
tree-pc-yuk.o: tree-pc-yuk.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
xmalloc.o: xmalloc.c xmalloc.h
sort.o: sort.c
main-direct.o: main-direct.c array.h xmalloc.h tools.h
txt2bin.o: txt2bin.c array.h xmalloc.h
makegrid.o: makegrid.c
