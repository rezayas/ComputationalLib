CFILES = lineartest.cpp polytest.cpp mintest.cpp lineartest2.cpp
OFILES = lineartest.o polytest.o mintest.o lineartest2.o
EFILES = lineartest polytest mintest lineartest2
all: $(EFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
lineartest.o: ../lib/linreg.hpp
lineartest2.o: ../lib/linreg.hpp
lineartest2: lineartest2.o ../lib/linreg.o
	c++ $^ -o lineartest2
lineartest: lineartest.o ../lib/linreg.o
	c++ $^ -o lineartest
polytest: polytest.o ../lib/linreg.o ../lib/polyreg.o ../lib/poly.o
	c++ $^ -o polytest
polytest.o: ../lib/polyreg.hpp
mintest.o: ../lib/poly.hpp
mintest: mintest.o ../lib/linreg.o ../lib/poly.o
	c++ $^ -o mintest
clean:
	rm -f $(OFILES) $(EFILES)
