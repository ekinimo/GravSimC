CC=gcc
CFLAGS= -Wall -Wextra -Og -ggdb3
DEPS= -lSDL2 -lm

sim: grav_sim.c
	$(CC) grav_sim.c $(CFLAGS) $(DEPS) -o simulation 
