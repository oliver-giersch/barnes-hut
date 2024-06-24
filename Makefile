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
SRC := src/main.c
INC := -I./include -I./src
LIB := -lpthread

# adjust SRC based on XWIN before!
OBJ := $(SRC:%.c=%.o)
