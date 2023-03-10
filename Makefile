all: typical-bezier-test

GEOM=../libgeom
INCLUDES=-I$(GEOM)
CXXFLAGS=-std=c++17 -Wall -pedantic -O3 -g -DNDEBUG $(INCLUDES)
LIBS=-L$(GEOM)/release -lgeom

typical-bezier-test: typical-bezier-test.o typical-bezier.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(INCLUDES) $(LIBS)
