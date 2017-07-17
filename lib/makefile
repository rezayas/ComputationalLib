OFILES = calc.o constr.o linreg.o no_lapack_workaround.o
HFILES = calc.hpp linreg.hpp no_lapack_workaround.hpp
CFILES = calc.cpp constr.cpp linreg.cpp no_lapack_workaround.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp $(if $(shell test -e %.hpp), %.hpp) no_lapack_workaround.hpp
	c++ -c -std=c++11 $<
clean:
	rm $(OFILES)
