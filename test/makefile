CFILES = lineartest.cpp polytest.cpp mintest.cpp lineartest2.cpp poly2test.cpp
OFILES = lineartest.o polytest.o mintest.o lineartest2.o poly2test.o
EFILES = lineartest polytest mintest linear2test poly2test
all: $(EFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp
	c++ -c -std=c++11 -I../lib $<
lineartest.o: ../lib/linreg.hpp
lineartest2.o: ../lib/linreg.hpp
linear2test: lineartest2.o ../lib/linreg.o
	c++ $^ -o linear2test
lineartest: lineartest.o ../lib/linreg.o
	c++ $^ -o lineartest
polytest: polytest.o ../lib/linreg.o ../lib/polyreg.o ../lib/poly.o
	c++ $^ -o polytest
polytest.o: ../lib/polyreg.hpp
poly2test: poly2test.o ../lib/linreg.o ../lib/polyreg.o ../lib/poly.o
	c++ $^ -o poly2test
poly2test.o: ../lib/polyreg.hpp
mintest.o: ../lib/poly.hpp
mintest: mintest.o ../lib/linreg.o ../lib/poly.o
	c++ $^ -o mintest
clean:
	rm -f $(OFILES) $(EFILES)
