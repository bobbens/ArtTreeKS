
# Options on what packages to use
MINPACK 	:= no  # deprecated fortran MINPACK package
CMINPACK := yes # modern C version of MINPACK (currently required)
NLOPT		:= no  # Proper mathematical optimization methods

# Objects to link
LIBNAME  := arttreeks
LIB		:= synthesis.so
LIB_OBJS := synthesis.o synth_lua.o solver_minpack.o solver_ga.o visualizer.o sminpack.o mapmask.o
LIB_OBJS += rand.o mem.o kin_misc.o
LIB_OBJS += cmaes.o solver_cmaes.o

# COMPILATION FLAGS
CFLAGS  := -g -O2 -fPIC
CFLAGS  += -W -Wall -Wextra -Wunused -Wshadow -Wmissing-prototypes -Winline -Wcast-align -Wmissing-declarations -Wredundant-decls -Wno-long-long -Wcast-align
#CFLAGS  += -DNDEBUG
CFLAGS  += $(shell pkg-config lua5.1 --cflags)
#CFLAGS  += $(shell pkg-config lua --cflags)
CFLAGS  += $(shell pkg-config sdl --cflags)
#CFLAGS += -DRAND_ACCURATE

# LINKER FLAGS
LDFLAGS := -ldq -lpthread -lcminpack
LDFLAGS += $(shell pkg-config lua5.1 --libs)
#LDFLAGS += $(shell pkg-config lua --libs)
LDFLAGS += $(shell pkg-config sdl --libs)
LDFLAGS += $(shell pkg-config gl --libs)
LDFLAGS += $(shell pkg-config glu --libs)
LDFLAGS += $(shell pkg-config libpng --libs)

# Fortran MINPACK conditional stuff
ifeq ($(MINPACK),yes)
LIB_OBJS += solver_fminpack.o
CFLAGS   += -DHAVE_MINPACK
LDFLAGS  += -lminpack
endif

# NLOPT conditional stuff
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
	$(CC) $(CFLAGS) -shared $(LIB_OBJS) -o $(LIB) $(LDFLAGS)

docs:
	doxygen
	$(MAKE) -C docs/latex
	cp docs/latex/refman.pdf $(LIBNAME).pdf




