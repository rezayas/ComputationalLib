OFILES = linreg.o polyreg.o poly.o qp.o multiquad.o csv_linreg.o
HFILES = linreg.hpp polyreg.hpp poly.hpp qp.hpp multiquad.hpp csv_linreg.hpp
CFILES = linreg.cpp polyreg.cpp poly.cpp qp.cpp multiquad.cpp csv_linreg.cpp
all: $(OFILES)
.PHONY: all clean
$(OFILES): %.o: %.cpp %.hpp
	c++ -c -std=c++11 $<
polyreg.hpp: linreg.hpp poly.hpp
	touch polyreg.hpp
poly.o: linreg.hpp
multiquad.hpp: polyreg.hpp
	touch multiquad.hpp
csv_linreg.o: linreg.hpp
clean:
	rm -f $(OFILES)
.PRECIOUS: $(HFILES)
