PROGNAME=endoev

IDIR =./src/include
ODIR=src/obj
SRCDIR=src

CC=g++

CFLAGS=-I$(IDIR) `pkg-config --cflags gsl`

LIBS=-lm `pkg-config --libs gsl` -lpng

_DEPS = dv_tools.h ca.h randomgen.h stringreps.h parameters.h mycol.h mypng.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o dv_tools.o ca.o stringreps.o parameters.o mycol.o mypng.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_test = test.o 
OBJ_test = $(patsubst %,$(ODIR)/%,$(_OBJ_test))

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

.PHONY: run

run:
	./$(PROGNAME)  

test: $(OBJ_test)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	echo $(OBJ_test) $(CFLAGS) $(LIBS)
	./test

