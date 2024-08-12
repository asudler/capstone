CXX = g++
CXXFLAGS = -std=c++17
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/src/rkf45/rkf45.cc \
	   ~/git/capstone/src/spline/spline.cc
OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

oscillator.dat : $(TARGET)
	./$< --test -b 0.25 -c 5.0 -yi 3.04159 -dyi 0 -ti 0 -tf 10 -dt 0.1 > $@

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf *.out *.o *.dat *.png

