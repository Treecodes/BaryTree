OFILES = main.o tools.o tree-cp.o tree-cp-yuk.o treedriver.o partition.o xmalloc.o
TARGET = treecp

OFILESD = main-direct.o tools.o xmalloc.o
TARGETD = directpot

CC = mpicc
CFLAGS = -Wall -Wextra  -O2 -std=c99
$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(TARGETD): $(OFILESD)
	$(CC) $(OFILESD) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(OFILESD) $(TARGETD)


main.o: main.c array.h xmalloc.h treedriver.h tools.h
tools.o: tools.c array.h xmalloc.h
treedriver.o: treedriver.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  tree-cp.h treedriver.h
partition.o: partition.c tools.h partition.h
tree-cp.o: tree-cp.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree-cp.h
xmalloc.o: xmalloc.c xmalloc.h
tree-cp-yuk.o: tree-cp-yuk.c array.h xmalloc.h globvars.h tnode.h tools.h \
		  partition.h tree-cp.h
main-direct.o: main-direct.c array.h xmalloc.h tools.h
