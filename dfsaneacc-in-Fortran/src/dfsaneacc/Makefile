all: qr.o dfsaneacc.o
	$(AR) rcs ./lib/libdfsaneacc.a $^
	$(shell if [ -e bmdfsaneacc.mod ]; then mv bmdfsaneacc.mod ./inc/; fi)
	$(shell if [ -e bmqr.mod ]; then mv bmqr.mod ./inc/; fi)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f dfsaneacc.o qr.o

distclean: clean
	rm -f ./inc/bmdfsaneacc.mod ./inc/bmqr.mod
	rm -f ./lib/libdfsaneacc.a
