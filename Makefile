EXE     := barnes-hut
CC      := /usr/bin/cc
CFLAGS  := --std=c11 -Wall -Wpedantic
LD      := /usr/bin/cc
LDFLAGS :=

COPTFLAGS := -O3 --fast-math
# COPTFLAGS := -O3  -fno-math-errno -fno-trapping-math

# structure:
# 	include/arena.h, inline function for malloc
# 	idea: phys.c instead of geom.h?
# 	src/geom.c <- include/geom.h (vec2 + inline fns, quadrant, particle_tree)
# 	src/main.c <- main, worker threads, synchronization, arena?
# 	src/options.c <- option parsing
SRC := src/main.c src/options.c src/phys.c
INC := -I./include
LIB := -lpthread -lmath

# adjust SRC based on XWIN before!
OBJ := $(SRC:%.c=%.o)
DEP := $(SRC:%.c=%.d)

all: $(EXE)

$(EXE): $(OBJ) Makefile
	$(LD) $(LDFLAGS) $< -o $@

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< $(INC) $(LIB) -o $@

clean:
	rm $(OBJ) $(DEP) $(EXE) 2> /dev/null || true

.PHONY: all clean
