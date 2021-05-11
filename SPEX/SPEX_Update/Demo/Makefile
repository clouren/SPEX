LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd
CS = $(LDLIBS)

all: test
	- ./test
test: test.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o test test.c $(CS)

