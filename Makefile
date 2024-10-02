BIN     := barnes-hut
CC      := /usr/bin/cc
CFLAGS  := --std=c11 -Werror -Wall -Wpedantic -MMD
LD      := /usr/bin/cc
LDFLAGS :=

# safer alternative: -O3 -fno-math-errno -fno-trapping-math
COPTFLAGS := -O3 -ffast-math

SRC := src/main.c src/options.c src/phys.c
INC := -I./include
LIB := -lpthread -lm

ifeq ($(RENDER),1)
	CFLAGS += -DRENDER
	SRC    += src/render.c
	LIB    += -lSDL2 -lGL -lGLU
endif

OBJ := $(SRC:%.c=%.o)
DEP := $(SRC:%.c=%.d)

ifneq (($BUILD),debug)
	CFLAGS  += $(COPTFLAGS) -flto
	LDFLAGS += -Wl,-flto
else
	CFLAGS += -g
endif

all: $(BIN)

compiledb: compile_commands.json

$(BIN): $(OBJ) Makefile
	$(LD) $(LDFLAGS) $(OBJ) $(LIB) -o $@

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

-include $(DEP)

clean:
	rm $(BIN) src/*.o src/*.d compile_commands.json 2> /dev/null || true

compile_commands.json:
	bear -- $(MAKE) RENDER=1 all

.PHONY: all compiledb clean
