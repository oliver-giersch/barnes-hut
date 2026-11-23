include options.mk

CFLAGS := --std=c23 -Werror -Wall -Wpedantic -MMD
CDEFS  := -D_DEFAULT_SOURCE

OPTFLAGS := -O3 -march=native -mtune=native -fno-math-errno -fno-trapping-math

SRC := src/main.c src/options.c src/physics.c src/thread.c
INC := -I./include
LIB := #-lpthread -lm? (-lhwloc)
OBJ := $(SRC:src/%.c=obj/%.o)
DEP := $(SRC:src/%.c=obj/%.d)

all: barnes-hut

barnes-hut: $(OBJ) Makefile
	$(LD) $(LDFLAGS) $(OBJ) $(LIB) -o $@

$(OBJ): obj/%.o: src/%.c | obj
	$(CC) $(CFLAGS) $(INC) $(CDEFS) -c $< -o $@

obj:
	mkdir -p $@

-include $(DEP)

.PHONY: all clean
