# Make file for Pairwise, Interval, Lkgen, Convert, fin and rhomap

CC = gcc
CPP = g++
CFLAGS = -O3 
#CFLAGS = -g -pedantic
LIB = -lm

OBJ = tools.o seqtools.o pairlk.o 

all: convert pairwise interval lkgen complete stat fin rhomap

convert: $(OBJ) convert.o tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS1) -o convert $(OBJ) convert.o $(LIB)

pairwise: $(OBJ) snp_sim.o pairwise.o tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS1) -o pairwise $(OBJ) snp_sim.o pairwise.o $(LIB)

interval: $(OBJ) pair_int.o tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS1) -o interval $(OBJ) pair_int.o $(LIB)

lkgen: $(OBJ) lkgen.o tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS) -o lkgen $(OBJ) lkgen.o $(LIB)

complete: $(OBJ) complete.o tools.h ldhat.h
	$(CC) $(CFLAGS) -o complete $(OBJ) complete.o $(LIB)

stat: tools.o stat.o tools.h
	$(CC) $(CFLAGS) -o stat tools.o stat.o $(LIB)
	
fin: fin.c gamma.c tools.c
	$(CC) $(CFLAGS) fin.c gamma.c tools.c $(LIB) -o fin
	
rhomap: rhomap.cpp compLK.cpp rhomap_tools.o 
	$(CPP) $(CFLAGS) rhomap.cpp compLK.cpp rhomap_tools.o $(LIB) -o rhomap

convert.o: $(OBJ) convert.c tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS) -c -o convert.o convert.c

pairwise.o:  $(OBJ) pairdip.c tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS) -c -o pairwise.o pairdip.c

pair_int.o: $(OBJ) pair_int.c tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS) -c -o pair_int.o pair_int.c

lkgen.o: $(OBJ) lkgen.c tools.h seqtools.h ldhat.h
	$(CC) $(CFLAGS) -c -o lkgen.o lkgen.c

complete.o: $(OBJ) complete.c tools.h ldhat.h
	$(CC) $(CFLAGS) -c -o complete.o complete.c

tools.o: tools.c tools.h
	$(CC) $(CFLAGS) -c -o tools.o tools.c

seqtools.o: seqtools.c seqtools.h ldhat.h
	$(CC) $(CFLAGS) -c -o seqtools.o seqtools.c 

pairlk.o: pairlk.c 
	$(CC) $(CFLAGS) -c -o pairlk.o pairlk.c

stat.o: stat.c
	$(CC) $(CFLAGS) -c -o stat.o stat.c

rhomap_tools.o: rhomap_tools.cpp
	$(CPP) $(CFLAGS) -c -o rhomap_tools.o rhomap_tools.cpp
		
clean:
	rm -rf *.o
	rm -rf *~
