CC      := /usr/bin/cc
CFLAGS  := --std=c11 -Wall -Wpedantic
LD      := /usr/bin/cc
LDFLAGS :=

COPTFLAGS := -O3 --fast-math
# COPTFLAGS := -O3  -fno-math-errno -fno-trapping-math

SRC := src/main.c
INC := -I./include -I./src
LIB := -lpthread

# adjust SRC based on XWIN before!
OBJ := $(SRC:%.c=%.o)
