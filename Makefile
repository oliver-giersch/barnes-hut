CFLAGS := --std=c23 -Werror -Wall -Wpedantic -MMD

SRC := src/main.c src/physics.c src/thread.c
INC := -I./include

OBJ := $(SRC:%.c=%.o)
DEP := $(SRC:%.c=%.d)

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) $(INC) $(CDEFS) -c $< -o $@
