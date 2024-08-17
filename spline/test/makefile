CXX = g++
CXXFLAGS = -std=c++17
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/misctools/misctools.cc \
       ~/git/capstone/spline/spline.cc
OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

out.cspline_im.png : in.cspline_im.gp out.cspline.csv in.fcn.csv
	gnuplot -e "outfile='$@' ; interpolation='$(word 2,$^)' ; \
	data='$(word 3,$^)'" $<

out.cspline_re.png : in.cspline_re.gp out.cspline.csv in.fcn.csv
	gnuplot -e "outfile='$@' ; interpolation='$(word 2,$^)' ; \
	data='$(word 3,$^)'" $<

out.cspline.csv : $(TARGET) in.fcn.csv
	./$< --test -ni 100 -fname $(word 2,$^) > $@

in.fcn.csv : $(TARGET)
	./$< --makedata -ni 12 -xmax 10 > $@

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean : 
	rm -rf *.out *.o *.dat *.png *.csv

