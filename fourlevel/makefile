CXX = g++
CXXFLAGS = -std=c++17 -g
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/fourlevel/fourlevel.cc \
	   ~/git/capstone/linalg/matrix.cc \
	   ~/git/capstone/misctools/misctools.cc

OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf out.* *.o *.png *.csv

