src     = addMsite.cpp
exes    = addMsite.exe
CC      = g++
LIBS    = -lgmx_reader -lxdrfile -lm

all: ${exes}

${exes}: ${src} addMsite.h
	$(CC) $(src) -o $(exes) $(LIBS) -std=c++11 -fmax-errors=10 -fopenmp -lpthread -O3

clean:
	rm addMsite.exe
