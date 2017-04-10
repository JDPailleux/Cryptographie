# MAKEFILE

DM.o : DM.c DM.h
	gcc -o DM.o -c DM.c -W -Wall -lm -lgmp

DM: DM.o
	gcc -o DM DM.o -lm -lgmp

all: DM

clean:
	rm -rf *.o
mrproper: clean
	rm -rf hello
run: 
	./DM

