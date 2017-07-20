CFILES = lineartest.cpp
OFILES = lineartest.o
all: lineartest
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
lineartest.o: ../lib/linreg.hpp
lineartest: lineartest.o ../lib/linreg.o
	c++ $^ -o lineartest
clean:
	-rm $(OFILES) lineartest
