CXX = g++
CXXFLAGS = -std=c++17 -g -fsanitize=address,undefined -fno-omit-frame-pointer
LDLIBS = -lstdc++

SRCS = main.cc ~/git/capstone/linalg/matrix/matrix.cc
OBJS = $(SRCS:.cc=.o)

TARGET = main.out
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) $^

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean :
	rm -rf *.out *.o *.dat *.png *.csv ../*.o

