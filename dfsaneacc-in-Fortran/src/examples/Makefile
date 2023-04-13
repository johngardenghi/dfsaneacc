all: dfsaneaccma.o
	$(FC) $(FFLAGS) -L$(DFSANEACCSRC) -o dfsaneaccma-example $^ -ldfsaneacc
	mkdir -p $(BIN)/examples
	mv dfsaneaccma-example $(BIN)/examples

dfsaneaccma.o: dfsaneaccma.f90
	$(FC) $(FFLAGS) -I$(DFSANEACCSRC) -c -o $@ $<

clean:
	rm -f dfsaneaccma.o

distclean: clean
	rm -f $(BIN)/examples/dfsaneaccma-example
