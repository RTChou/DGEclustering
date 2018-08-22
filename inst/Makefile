CC = gcc

LDLIBS = -lm

CFLAGS = -ggdb

PROGRAM = genclust

SOURCES.c = genclust.c genclustlb.c

OBJ=$(SOURCES.c:%.c=%.o)

$(PROGRAM): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDLIBS) -o $@ 

%.o:%.c
	$(CC) $(CFLAGS) -c $<
