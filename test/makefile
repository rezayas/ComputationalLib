CFILES = testa.cpp testb.cpp
OFILES = testa.o testb.o
all: testa testb
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
testa: testa.o
	c++ testa.o ../lib/calc.o -larmadillo -o testa
testb: testb.o
	c++ testb.o ../lib/linreg.o -larmadillo -o testb
clean:
	rm $(OFILES) testa testb
