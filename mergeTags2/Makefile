CC=g++
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=kstring.h alglib/*.h
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz -lm

mergeTags2: mergeTags2.cpp $(HEADERS) $(OBJECTS) dna.h
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)
