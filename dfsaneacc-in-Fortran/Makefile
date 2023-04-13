AR := ar
FC := gfortran

FFLAGS := -O3

DFSANEACC     := $(CURDIR)
BIN          := $(DFSANEACC)/bin
SRC          := $(DFSANEACC)/src
DFSANEACCSRC := $(SRC)/dfsaneacc
EXAMPLESSRC  := $(SRC)/examples
INTERSRC     := $(SRC)/interfaces

export

all:
	make -C $(DFSANEACCSRC) all

cutest: all
	make -C $(INTERSRC)/cutest all

cutest-distclean:
	make -C $(INTERSRC)/cutest distclean

example: all
	make -C $(EXAMPLESSRC) all

clean:
	make -C $(DFSANEACCSRC) clean
	make -C $(INTERSRC)/cutest clean
	make -C $(EXAMPLESSRC) clean

distclean: clean
	make -C $(DFSANEACCSRC) distclean
	make -C $(INTERSRC)/cutest distclean
	make -C $(EXAMPLESSRC) distclean
