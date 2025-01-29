CXX = g++
CXXFLAGS = -std=c++17 -g -fsanitize=address
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/fourlevel/fourlevel.cc \
	   ~/git/capstone/linalg/matrix/matrix.cc \
	   ~/git/capstone/misctools/misctools.cc \
	   ~/git/capstone/fourlevel/spatial.cc \
	   ~/git/capstone/spline/spline.cc

OBJS = $(SRCS:.cc=.o)

TARGET = main.out

out.rho.dat : $(TARGET) in.parameters.dat
	./$< -i in.parameters.dat > $@

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf out.* *.o *.png *.csv log.* *.log ../*.o

