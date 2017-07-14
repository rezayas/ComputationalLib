OFILES = calc.o constr.o linreg.o
HFILES = calc.hpp linreg.hpp
CFILES = calc.cpp constr.cpp linreg.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp $(if $(shell test -e %.hpp), %.hpp)
	c++ -c -std=c++11 $<
clean:
	rm $(OFILES)
