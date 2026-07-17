CXX ?= g++
PYTHON ?= python3

CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
OPENMP_FLAGS ?= -fopenmp

ifeq ($(OS),Windows_NT)
EXEEXT := .exe
else
EXEEXT :=
endif

TARGET := ngia3$(EXEEXT)
SOURCES := main.cpp makedb.cpp
HEADERS := makedb.hpp timer.hpp

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(SOURCES) -o $@

test: $(TARGET)
	cd .. && NGIA3_BINARY="$(CURDIR)/$(TARGET)" $(PYTHON) -m unittest discover -s tests -v

clean:
	$(RM) $(TARGET)
