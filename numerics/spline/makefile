CXX = g++
CXXFLAGS = -std=c++23
LDLIBS = -lstdc++

SRCS = main.cc misctools.cc spline.cc
OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

cspline.png : cspline.gp cspline.dat fcn.dat
	gnuplot -e "outfile='$@' ; interpolation='$(word 2,$^)' ; \
	data='$(word 3,$^)'" $<

cspline.dat : $(TARGET) fcn.dat
	./$< --test -ni 100 -fname $(word 2,$^) > $@

fcn.dat : $(TARGET)
	./$< --makedata -ni 12 -xmax 10 > $@

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean : 
	rm -rf *.out *.o *.dat *.png

