


MINPACK 	:= no
CMINPACK := yes
NLOPT		:= no


LIB		:= synthesis.so
LIB_OBJS := synthesis.o synth_lua.o solver_minpack.o solver_ga.o visualizer.o sminpack.o mapmask.o
LIB_OBJS += rand.o mem.o kin_misc.o
LIB_OBJS += cmaes.o solver_cmaes.o

CFLAGS  := -g -fPIC
CFLAGS  += -W -Wall -Wextra -Wunused -Wshadow -Wmissing-prototypes -Winline -Wcast-align -Wmissing-declarations -Wredundant-decls -Wno-long-long -Wcast-align
#CFLAGS  += -DNDEBUG
CFLAGS  += $(shell pkg-config lua5.1 --cflags)
CFLAGS  += $(shell pkg-config sdl --cflags)
LDFLAGS += $(shell pkg-config libpng --cflags)
#CFLAGS  += $(shell pkg-config lua --cflags)
#CFLAGS += -DRAND_ACCURATE
LDFLAGS := -ldq -lpthread -lcminpack
LDFLAGS += $(shell pkg-config lua5.1 --libs)
LDFLAGS += $(shell pkg-config sdl --libs)
LDFLAGS += $(shell pkg-config gl --libs)
LDFLAGS += $(shell pkg-config glu --libs)
LDFLAGS += $(shell pkg-config libpng --libs)
#LDFLAGS += $(shell pkg-config lua --libs)

ifeq ($(MINPACK),yes)
LIB_OBJS += solver_fminpack.o
CFLAGS   += -DHAVE_MINPACK
LDFLAGS  += -lminpack
endif

ifeq ($(NLOPT),yes)
LIB_OBJS += solver_nlopt.o
CFLAGS   += -DHAVE_NLOPT
LDFLAGS  += -lnlopt
endif

OBJS 		:= $(LIB_OBJS)

.PHONY: all clean


all: $(APP) $(LIB)

clean:
	$(RM) $(OBJS) $(APP)

$(LIB): $(LIB_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared $(LIB_OBJS) -o $(LIB)




