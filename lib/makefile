OFILES = linreg.o polyreg.o poly.o
HFILES = linreg.hpp polyreg.hpp poly.hpp
CFILES = linreg.cpp polyreg.cpp poly.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp %.hpp
	c++ -c -std=c++11 $<
polyreg.hpp: linreg.hpp poly.hpp
	touch polyreg.hpp
clean:
	-rm $(OFILES)
