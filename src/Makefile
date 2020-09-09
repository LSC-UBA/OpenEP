CC=g++

# Standar compilation
CFLAGS_STD=-c -Wall -Wextra 
LDFLAGS_STD=

# With O3 optimizations
CFLAGS_OPT=-c -O3 -Wall -Wextra 
LDFLAGS_OPT=

# With Ofast optimizations
CFLAGS_FAST=-c -Wall -Wextra  -Ofast
LDFLAGS_FAST=

# Debug compilation
CFLAGS_DBG=-c -g -Wall -Wextra  
LDFLAGS_DBG=

# With OpenMP
CFLAGS_OMP=-c -fopenmp -Wall -Wextra 
LDFLAGS_OMP= -fopenmp 

# With Ofast optimizations and OpenMP
CFLAGS_FAST_OMP=-c -Ofast -fopenmp -Wall -Wextra 
LDFLAGS_FAST_OMP= -fopenmp 

# With 03 optimizations and OpenMP
CFLAGS_O3_OMP=-c -O3 -fopenmp -Wall -Wextra 
LDFLAGS_O3_OMP= -fopenmp 

COMP_TYPE=FAST_OMP

CFLAGS  = ${CFLAGS_${COMP_TYPE}} --std=c++17
LDFLAGS = ${LDFLAGS_${COMP_TYPE}}

SOURCES=mesh.cpp scalar_field.cpp vector_field.cpp electrics_calc.cpp temp_calc.cpp save.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=simulate

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ 
	
clean:
	rm $(EXECUTABLE) -f
	rm *.o *.gch -f
