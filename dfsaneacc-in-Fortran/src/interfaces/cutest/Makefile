SIF := ELFUN EXTER GROUP RANGE
SIFPROB := $(shell if [ -e $(MASTSIF)/$(PROBNAME).SIF ]; then echo true; fi)

ifeq ($(SIFPROB), true)

sifdec:
	$(SIFDECODE)/bin/sifdecoder -A $(MYARCH) $(PROBNAME).SIF
	$(FC) $(FFLAGS) -c $(addsuffix .f, $(SIF))

else

sifdec:
	$(error $(PROBNAME).SIF does not exist in MASTSIF=$(MASTSIF) directory.)

endif

all: sifdec link

link: dfsaneaccma.o $(addsuffix .o, $(SIF))
	$(FC) $(FFLAGS) -L$(DFSANEACCSRC) -L$(CUTEST)/objects/$(MYARCH)/double -o dfsaneaccma-forcutest $^ -lcutest -ldfsaneacc
	mkdir -p $(BIN)/cutest
	mv OUTSDIF.d $(BIN)/cutest
	mv dfsaneaccma-forcutest $(BIN)/cutest

dfsaneaccma.o: dfsaneaccma.f90
	$(FC) $(FFLAGS) -I$(DFSANEACCSRC) -c -o $@ $<

clean:
	rm -f dfsaneaccma.o
	rm -f $(addsuffix .o, $(SIF))
	rm -f $(addsuffix .f, $(SIF))
	rm -f AUTOMAT.d

distclean: clean
	rm -f $(BIN)/cutest/dfsaneaccma-forcutest
	rm -f $(BIN)/cutest/OUTSDIF.d
