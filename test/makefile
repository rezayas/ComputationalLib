CFILES = testb.cpp
OFILES = testb.o
all: testb
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
testb: testb.o
	c++ testb.o ../lib/linreg.o -larmadillo -o testb
clean:
	-rm $(OFILES) testb
