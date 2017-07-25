OFILES = linreg.o polyreg.o poly.o qp.o
HFILES = linreg.hpp polyreg.hpp poly.hpp qp.hpp
CFILES = linreg.cpp polyreg.cpp poly.cpp qp.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp %.hpp
	c++ -c -std=c++11 $<
polyreg.hpp: linreg.hpp poly.hpp
	touch polyreg.hpp
poly.o: linreg.hpp
clean:
	rm -f $(OFILES)
.PRECIOUS: polyreg.hpp
