#############################################
# Makefile to generate the python module
#############################################

CC=g++
CFLAGS=-shared -fPIC -O3
INCLUDES=$(shell python-config --includes) $(shell lhapdf-config --cflags) $(shell applgrid-config --cxxflags)
LIBS=$(shell python-config --libs) $(shell lhapdf-config --libs) $(shell applgrid-config --ldflags)

all: release

release: applwrap.so

applwrap.so: applwrap.cc
	$(CC) $(CFLAGS) $(INCLUDES) $< $(LIBS) -o $@

clean:
	rm -rf applwrap.so *~