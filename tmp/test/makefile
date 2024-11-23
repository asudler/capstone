CXX = g++
CXXFLAGS = -std=c++17 -g
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/tmp/fourlevel.cc \
	   ~/git/capstone/linalg/matrix.cc \
	   ~/git/capstone/misctools/misctools.cc \
	   ~/git/capstone/spline/spline.cc
	   #~/git/capstone/rkf45/rkf45.cc \

OBJS = $(SRCS:.cc=.o)

TARGET = main.out

out.rho.dat : $(TARGET) in.inputfile.dat
	./$< -i in.inputfile.dat > $@

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf out.* *.o *.png *.csv log.* *.log ../*.o

