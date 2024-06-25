EXE     := barnes-hut
CC      := /usr/bin/cc
CFLAGS  := --std=c11 -Werror -Wall -Wpedantic
LD      := /usr/bin/cc
LDFLAGS := -Wl,-lm

# safer alternative: -O3 -fno-math-errno -fno-trapping-math
COPTFLAGS := -O3 --fast-math

SRC := src/main.c src/options.c src/phys.c
INC := -I./include
LIB := -lpthread -lm

# adjust SRC based on XWIN before!
OBJ := $(SRC:%.c=%.o)
DEP := $(SRC:%.c=%.d)

ifeq ($(BUILD),optimize)
	CFLAGS  += $(COPTFLAGS)
else ifeq ($(BUILD),optimize-lto)
	CFLAGS  += $(COPTFLAGS) -flto
	LDFLAGS += -Wl,-flto
endif

all: $(EXE)

$(EXE): $(OBJ) Makefile
	$(LD) $(LDFLAGS) $(OBJ) -o $@

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INC) $(LIB)

clean:
	rm $(OBJ) $(DEP) $(EXE) 2> /dev/null || true

.PHONY: all clean
