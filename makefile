# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -std=c++11 -fopenmp

# the build target executable:
TARGET = main

all: $(TARGET)

# Use -DDEBUG=1, 2 or 3. The high DEBUG the more verbose the program.

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -DDEBUG=2 -o $(TARGET) $(TARGET).cpp Array_3D.cpp Materials.cpp

clean:
	$(RM) $(TARGET)
