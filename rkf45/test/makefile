CXX = g++
CXXFLAGS = -std=c++17 -g
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/spline/spline.cc
OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

out.oscillator.png : in.oscillator.gp out.oscillator.csv
	gnuplot -e "outfile='$@' ; data='$(word 2,$^)'" $<

out.oscillator.csv : $(TARGET)
	./$< --test -b 0.25 -c 5.0 -yi 3.04159 -dyi 0 -ti 0 -tf 10 -dt 0.1 > $@

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf *.out *.o *.dat *.png *.csv

