
CXX := clang++
CXXFLAGS := -O3 -Wall -Wextra -std=c++11 -march=native

all:
	$(CXX) $(CXXFLAGS) main.cpp -I. -o DelaunayTri
