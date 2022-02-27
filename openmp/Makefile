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
CENTROIDS   =   50
ELEMENTS    =   50000

KMEANS_INPUT_FILE = input.txt
KMEANS_OUTPUT_FILE = output.txt

# all : $(EXE)
# 	@echo "Building..."

# $(EXE) : $(OBJECTS)
# 	$(LD) $^ $(CFLAGS) $(LDFLAGS) -o $@ -fopenmp

# %.o : %.c
# 	$(CC) $(CFLAGS) -c $< -o $@

# %.o : %.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@


build_open_mp :
	gcc -o kmeans kmeans.c -fopenmp

build_sequential :
	gcc -o kmeans kmeans.c

clean :
	rm -fr $(EXE) $(OBJECTS)

gen :
	python3 $(PYTHONFILE) $(CENTROIDS) $(ELEMENTS) > $(KMEANS_INPUT_FILE)

run : 
	./kmeans < $(KMEANS_INPUT_FILE) > $(KMEANS_OUTPUT_FILE)

all : $(clean) $(build_open_mp) $(gen) $(run)
	@echo "Running"
	
# run : $(clean) $(EXE) 
# 	./$(EXE) < $(KMEANS_INPUT_FILE) > $(KMEANS_OUTPUT_FILE)


