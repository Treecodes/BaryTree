OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o tree-pc.o tree-pc-yuk.o \
	 treedriver.o partition.o xmalloc.o sort.o
TARGET = tree

OFILESD = main-grid.o tools.o tree-cp.o tree-cp-grid.o \
	  treedriver-grid.o partition.o xmalloc.o sort.o
TARGETD = tree-grid

OFILESE = main-grid-bdry.o tools.o tree-cp.o tree-cp-grid-bdry.o \
      treedriver-grid-bdry.o partition.o xmalloc.o sort.o
TARGETE = tree-grid-bdry

OFILESA = main-direct.o tools.o xmalloc.o
TARGETA = direct

OFILESB = txt2bin.o sort.o xmalloc.o
TARGETB = txt2bin

OFILESC = makegrid.o
TARGETC = makegrid

CC = mpicc
CFLAGS = -Wall -O3
$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(OFILES) -o $@

$(TARGETA): $(OFILESA)
	$(CC) $(CFLAGS) $(OFILESA) -o $@

$(TARGETB): $(OFILESB)
	$(CC) $(CFLAGS) $(OFILESB) -o $@

$(TARGETC): $(OFILESC)
	$(CC) $(CFLAGS) $(OFILESC) -o $@

$(TARGETD): $(OFILESD)
	$(CC) $(CFLAGS) $(OFILESD) -o $@

$(TARGETE): $(OFILESE)
	$(CC) $(CFLAGS) $(OFILESE) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) $(OFILESE) $(TARGETE)
	rm -f $(OFILESB) $(TARGETB) $(TARGETC) $(OFILESC) $(TARGETD) $(OFILESD)$


main.o: main.c array.h xmalloc.h treedriver.h tools.h
main-grid.o: main-grid.c array.h xmalloc.h treedriver.h tools.h
main-grid-bdry.o: main-grid-bdry.c array.h xmalloc.h treedriver.h tools.h
tools.o: tools.c array.h xmalloc.h
treedriver.o: treedriver.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  tree.h treedriver.h
treedriver-grid.o: treedriver-grid.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  tree.h treedriver.h
treedriver-grid-bdry.o: treedriver-grid-bdry.c array.h xmalloc.h globvars.h tnode.h tools.h \
          tree.h treedriver.h
partition.o: partition.c tools.h partition.h
tree-cp.o: tree-cp.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
tree-cp-grid.o: tree-cp-grid.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree.h
tree-cp-grid-bdry.o: tree-cp-grid-bdry.c array.h xmalloc.h globvars.h tnode.h tools.h \
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
