CPP=g++

#PROFILE=-pg
#DEBUG=-DDEBUG
ifndef DEBUG
OPT=-Ofast -pipe -m64 -fwhole-program -march=native -mtune=native # -flto ("do not use with -fwhole-program, according to gcc manual")
endif

GDB=-ggdb3
LIBS=-lm -lsdpa -llapack -lblas -lgfortran -lquadmath -lpthread -ldmumps_seq
STATIC=-static-libgcc -static-libstdc++
WARNS=-Wextra -Wall -Wno-sign-compare -Wshadow -Wstrict-aliasing=1 -Werror -Wno-unused-result -pedantic-errors

CFLAGS=-ansi -std=c++11 -fabi-version=6 $(WARNS) $(OPT) $(STATIC) $(LIBS) $(GDB) $(DEBUG) $(PROFILE)

GRASS_DEPS=cutnorm.h graph.h grass.h summ.h util.h Makefile
SUMM_DEPS=cutnorm.h graph.h szem_clust.h szem_fk.h summ.h util.h sparse.h fenwick.h Makefile 

.PHONY: all clean

all: gen grass queries random sample_graph summ summ_to_dense_matrix

clean:
	-/bin/rm -rf gen data/gen grass queries random sample_graph summ summ_to_dense_matrix *.dSYM

gen: data/gen.cpp util.h Makefile
	$(CPP) data/$@.cpp -std=c++11 $(OPT) $(WARNS) $(STATIC) $(GDB) -o data/$@
	ln -fs data/gen gen

grass: grass.cpp $(GRASS_DEPS)
	$(CPP) $@.cpp $(CFLAGS) -o $@ 

queries: queries.cpp $(GRASS_DEPS)
	$(CPP) $@.cpp $(CFLAGS) -o $@

random: random.cpp $(GRASS_DEPS)
	$(CPP) $@.cpp $(CFLAGS) -o $@ 

sample_graph: sample_graph.cpp graph.h util.h Makefile
	$(CPP) $@.cpp $(CFLAGS) -o $@

summ: summ.cpp $(SUMM_DEPS)
	$(CPP) $@.cpp $(CFLAGS) -o $@ 

summ_to_dense_matrix: summ_to_dense_matrix.cpp $(SUMM_DEPS)
	$(CPP) $@.cpp $(CFLAGS) -o $@ 

# g++ summ.cpp -std=c++11 -lm -lsdpa -llapack -lblas -lgfortran -lquadmath -lpthread -ldmumps_seq -o summ
