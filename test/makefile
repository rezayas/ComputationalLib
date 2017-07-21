CFILES = lineartest.cpp polytest.cpp
OFILES = lineartest.o polytest.o
EFILES = lineartest polytest
all: $(EFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
lineartest.o: ../lib/linreg.hpp
lineartest: lineartest.o ../lib/linreg.o
	c++ $^ -o lineartest
polytest: polytest.o ../lib/linreg.o ../lib/polyreg.o ../lib/poly.o
	c++ $^ -o polytest
polytest.o: ../lib/polyreg.hpp
clean:
	rm -f $(OFILES) $(EFILES)
