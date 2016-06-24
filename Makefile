VERSION = 1.3.1
DATE = Jan 11, 2005

VER = "\"$(VERSION); $(DATE)\""
CMP = gcc -c
CMP_FLAGS =
LD = gcc 
LD_FLAGS = -lm
UNIFORM_RANDOM = drand48
RANDOM_SEED = srand48

snphap: snphap.o snphapfun.o cline.o 
	$(LD) $(LD_FLAGS)  -o snphap snphap.o snphapfun.o cline.o 

snphap.o : snphap.c snphap.h cline.h 
	$(CMP) $(CMP_FLAGS)  -DSEED=$(RANDOM_SEED) -DVERSION=$(VER) snphap.c

snphapfun.o : snphapfun.c  snphap.h
	$(CMP) $(CMP_FLAGS) -DURAN=$(UNIFORM_RANDOM) snphapfun.c

cline.o : cline.c cline.h
	$(CMP) $(CMP_FLAGS) cline.c

clean :
	rm -f *.o *~ snphap-$(VERSION).tar.gz snphap-$(VERSION).zip

distribution:
	mkdir snphap-$(VERSION)
	cp *.c *.h *.doc test.dat test.nam Makefile  snphap-$(VERSION)/
	rm -f snphap-$(VERSION).tar.gz snphap-$(VERSION).zip
	tar cvf snphap-$(VERSION).tar snphap-$(VERSION)/*
	gzip snphap-$(VERSION).tar
	zip snphap-$(VERSION).zip snphap-$(VERSION)/*
	rm -f snphap-$(VERSION)/*
	rmdir snphap-$(VERSION)





