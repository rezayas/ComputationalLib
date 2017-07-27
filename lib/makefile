OFILES = linreg.o polyreg.o poly.o qp.o multiquad.o
HFILES = linreg.hpp polyreg.hpp poly.hpp qp.hpp multiquad.hpp
CFILES = linreg.cpp polyreg.cpp poly.cpp qp.cpp multiquad.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp %.hpp
	c++ -c -std=c++11 $<
polyreg.hpp: linreg.hpp poly.hpp
	touch polyreg.hpp
poly.o: linreg.hpp
multiquad.hpp: polyreg.hpp
	touch multiquad.hpp
clean:
	rm -f $(OFILES)
.PRECIOUS: polyreg.hpp
