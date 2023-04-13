all: qr.o dfsaneacc.o
	$(AR) rcs libdfsaneacc.a $^

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f dfsaneacc.o qr.o

distclean: clean
	rm -f bmdfsaneacc.mod bmqr.mod
	rm -f libdfsaneacc.a
