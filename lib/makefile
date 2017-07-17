OFILES = linreg.o
HFILES = linreg.hpp
CFILES = linreg.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp $(if $(shell test -e %.hpp), %.hpp)
	c++ -c -std=c++11 $<
clean:
	-rm $(OFILES)
