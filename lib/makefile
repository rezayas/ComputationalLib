OFILES = linreg.o polyreg.o
HFILES = linreg.hpp polyreg.hpp
CFILES = linreg.cpp polyreg.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp $(wildcard %.hpp)
	c++ -c -std=c++11 $<
polyreg.hpp: linreg.hpp
	touch polyreg.hpp
clean:
	-rm $(OFILES)
