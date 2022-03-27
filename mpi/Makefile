NAME		=	kmeans
CC          =   gcc
LD			=	$(CC)
DEBUG		+=	#-DDEBUG
CFLAGS		+=	$(DEBUG)
CXXFLAGS	+=	$(CFLAGS)
LDFLAGS		+=

EXE			=	$(NAME)
CFILES		=	$(wildcard *.c)
CXXFILES	=	$(wildcard *.cpp)
OBJECTS		=	$(CFILES:.c=.o) $(CXXFILES:.cpp=.o)

PYTHONFILE  =   ./geninput.py
CENTROIDS   =   20
ELEMENTS    =   1000
NUM_THREADS =   2

KMEANS_INPUT_FILE = input.txt
KMEANS_OUTPUT_FILE_SEQ = output-times-seq.txt
KMEANS_OUTPUT_FILE_PAR = output-times-par.txt
KMEANS_TIME_MEASURE = time-measures.txt

# all : $(EXE)
# 	@echo "Building..."

# $(EXE) : $(OBJECTS)
# 	$(LD) $^ $(CFLAGS) $(LDFLAGS) -o $@ -fopenmp

# %.o : %.c
# 	$(CC) $(CFLAGS) -c $< -o $@

# %.o : %.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

gen :
	python3 $(PYTHONFILE) $(CENTROIDS) $(ELEMENTS) > $(KMEANS_INPUT_FILE)

run:
	# @rm ./kmeans
	@echo "Building..."
	@mpicc kmeans.c -o kmeans
	@echo "Running..."
	@mpirun -np 4 ./kmeans < $(KMEANS_INPUT_FILE)

clean :
	rm -fr $(EXE) $(OBJECTS)
	
# run : $(clean) $(EXE) 
# 	./$(EXE) < $(KMEANS_INPUT_FILE) > $(KMEANS_OUTPUT_FILE)


