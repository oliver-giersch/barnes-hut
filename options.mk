ifeq ($(OPTION_RENDER),1)
    OPTION_CDEFS += -DOPTION_RENDER
    OPTION_SRC   += src/render.c
    OPTION_LIB   += -lSDL3 -lGL
endif

ifeq ($(OPTION_NUMA),1)
    OPTION_CDEFS += -DOPTION_NUMA
    OPTION_SRC   += src/numa.c
    OPTION_LIB   += -lhwloc
endif
