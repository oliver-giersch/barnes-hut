BIN     := barnes-hut
CC      := /usr/bin/cc
CFLAGS  := --std=c11 -Werror -Wall -Wpedantic -MMD
LD      := /usr/bin/cc
LDFLAGS :=

# safer alternative: -O3 -fno-math-errno -fno-trapping-math
COPTFLAGS := -O3 --fast-math

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

ifeq ($(BUILD),optimize)
	CFLAGS  += $(COPTFLAGS)
else ifeq ($(BUILD),optimize-lto)
	CFLAGS  += $(COPTFLAGS) -flto
	LDFLAGS += -Wl,-flto
else
	CFLAGS  += -g
endif

all: $(BIN)

$(BIN): $(OBJ) Makefile
	$(LD) $(LDFLAGS) $(OBJ) $(LIB) -o $@

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

-include $(DEP)

clean:
	rm $(BIN) src/*.o src/*.d 2> /dev/null || true

.PHONY: all clean
