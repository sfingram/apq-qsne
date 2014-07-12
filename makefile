P=testapq 
OBJECTS=clustertree.o apq.o vec.o disjointset.o pqueue.o quadtree.o
#OBJECTS=clustertree.o apq.o vec.o disjointset.o pqueue.o
CFLAGS = `pkg-config --cflags glib-2.0` -g -Wall -O3
CXXFLAGS =`pkg-config --cflags glib-2.0` -g -Wall -O3 
LDLIBS=`pkg-config --libs glib-2.0` 
CC=g++
CXX=g++

$(P): $(OBJECTS)

prec-recall: prec-recall.o

tsne: quadtree.o vec.o clustertree.o disjointset.o

csv2dmat: vec.o csv.o

layoutquality: vec.o

clean:
	rm -rf *.o