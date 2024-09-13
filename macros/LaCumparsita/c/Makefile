GCC = gcc -ggdb 

OBJECTS = america.o \
	contor.o \
	hardspheres.o  \
	mountain2.o \
	piecefit.o \
	ellipsoid.o \
	bratu3Db.o \
	kissing2.o \
	packccmn.o \
	simfock2.o \
	cache.o \
	genpack-cc-mina.o \
	kissing.o \
	packccmn-feas.o \
	simfock.o \
	condor.o \
	genpack-csq-mina.o \
	location.o \
	packcrmn-feas.o \
	contor2.o \
	hardcube.o \
	mountain1.o \
	pedazos4.o \
	drand.o


all: $(OBJECTS) sci_lacumparsita.so

america.o: america.c
	$(GCC) -c america.c -o america.o
contor.o: contor.c
	$(GCC) -c contor.c -o contor.o
hardspheres.o: hardspheres.c
	$(GCC) -c hardspheres.c -o hardspheres.o
mountain2.o: mountain2.c
	$(GCC) -c mountain2.c -o mountain2.o
piecefit.o: piecefit.c
	$(GCC) -c piecefit.c -o piecefit.o
bratu3Db.o: bratu3Db.c
	$(GCC) -c bratu3Db.c -o bratu3Db.o
ellipsoid.o: ellipsoid.c
	$(GCC) -c ellipsoid.c -o ellipsoid.o
kissing2.o: kissing2.c
	$(GCC) -c kissing2.c -o kissing2.o
packccmn.o: packccmn.c
	$(GCC) -c packccmn.c -o packccmn.o
simfock2.o: simfock2.c
	$(GCC) -c simfock2.c -o simfock2.o
cache.o: cache.c
	$(GCC) -c cache.c -o cache.o
genpack-cc-mina.o: genpack-cc-mina.c
	$(GCC) -c genpack-cc-mina.c -o genpack-cc-mina.o
kissing.o: kissing.c
	$(GCC) -c kissing.c -o kissing.o
packccmn-feas.o: packccmn-feas.c
	$(GCC) -c packccmn-feas.c -o packccmn-feas.o
simfock.o: simfock.c
	$(GCC) -c simfock.c -o simfock.o
condor.o: condor.c
	$(GCC) -c condor.c -o condor.o
genpack-csq-mina.o: genpack-csq-mina.c
	$(GCC) -c genpack-csq-mina.c -o genpack-csq-mina.o
location.o: location.c
	$(GCC) -c location.c -o location.o
packcrmn-feas.o: packcrmn-feas.c
	$(GCC) -c packcrmn-feas.c -o packcrmn-feas.o
contor2.o: contor2.c
	$(GCC) -c contor2.c -o contor2.o
hardcube.o: hardcube.c
	$(GCC) -c hardcube.c -o hardcube.o
mountain1.o: mountain1.c
	$(GCC) -c mountain1.c -o mountain1.o
pedazos4.o: pedazos4.c
	$(GCC) -c pedazos4.c -o pedazos4.o
drand.o: drand.c
	$(GCC) -c drand.c -o drand.o

sci_lacumparsita.so: sci_lacumparsita.c
	$(GCC) -shared -o sci_lacumparsita.so sci_lacumparsita.c *.o -I/local/stow/scilab-4.1.1/lib/scilab-4.1.1/routines /usr/lib/libf2c.a

clean:
	rm *.o
