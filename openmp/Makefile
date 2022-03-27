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
CENTROIDS   =   100
ELEMENTS    =   100000
NUM_THREADS =   40

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



build_all:	
	sed -i 's/INPUT_THREADS/$(NUM_THREADS)/g' kmeans_atomic.c
	sed -i 's/INPUT_THREADS/$(NUM_THREADS)/g' kmeans_critical.c
	sed -i 's/INPUT_THREADS/$(NUM_THREADS)/g' kmeans_reduction.c
	gcc -o ./parallel/kmeans_atomic kmeans_atomic.c -fopenmp
	gcc -o ./parallel/kmeans_critical kmeans_critical.c -fopenmp
	gcc -o ./parallel/kmeans_reduction kmeans_reduction.c -fopenmp
	gcc -o ./seq/kmeans kmeans_atomic.c
	sed -i 's/NUM_THREADS $(NUM_THREADS)/NUM_THREADS INPUT_THREADS/g' kmeans_atomic.c
	sed -i 's/NUM_THREADS $(NUM_THREADS)/NUM_THREADS INPUT_THREADS/g' kmeans_critical.c
	sed -i 's/NUM_THREADS $(NUM_THREADS)/NUM_THREADS INPUT_THREADS/g' kmeans_reduction.c

clean_all:
	rm -rf ./parallel/kmean**
	rm -rf ./seq/kmeans
	rm -rf ./parallel/kmeans


build_open_mp :
	gcc -o ./parallel/kmeans kmeans.c -fopenmp

build_sequential :
	gcc -o ./seq/kmeans kmeans.c

clean :
	rm -fr $(EXE) $(OBJECTS)

gen :
	python3 $(PYTHONFILE) $(CENTROIDS) $(ELEMENTS) > $(KMEANS_INPUT_FILE)

run :
	python3 $(PYTHONFILE) $(CENTROIDS) $(ELEMENTS) > $(KMEANS_INPUT_FILE)
	sed -i 's/INPUT_THREADS/$(NUM_THREADS)/g' kmeans.c

	rm -fr kmeans kmeans.o

	gcc -o kmeans -fopenmp kmeans.c
	for number in 0 1 2; do \
		./kmeans < $(KMEANS_INPUT_FILE) | grep All >> $(KMEANS_OUTPUT_FILE_PAR); \
	done

	rm -fr kmeans kmeans.o

	gcc -o kmeans kmeans.c
	for number in 0 1 2; do \
		./kmeans < $(KMEANS_INPUT_FILE) | grep All >> $(KMEANS_OUTPUT_FILE_SEQ); \
	done

	rm -fr kmeans kmeans.o

	echo "Sequential Num threads: 1, num_elements: $(ELEMENTS), time measured: " >> $(KMEANS_TIME_MEASURE)
	python3 calculate_median.py -f $(KMEANS_OUTPUT_FILE_SEQ) >> $(KMEANS_TIME_MEASURE)
	echo "\n" >> $(KMEANS_TIME_MEASURE)

	echo "Parallel Num threads: $(NUM_THREADS), num_elements: $(ELEMENTS), time measured: " >> $(KMEANS_TIME_MEASURE)
	python3 calculate_median.py -f $(KMEANS_OUTPUT_FILE_PAR) >> $(KMEANS_TIME_MEASURE)
	echo "\n" >> $(KMEANS_TIME_MEASURE)

	rm $(KMEANS_OUTPUT_FILE_SEQ)
	rm $(KMEANS_OUTPUT_FILE_PAR)
	sed -i 's/NUM_THREADS $(NUM_THREADS)/NUM_THREADS INPUT_THREADS/g' kmeans.c

	cat $(KMEANS_TIME_MEASURE)

test:
	./kmeans < $(KMEANS_INPUT_FILE)

time :
	python3 $(PYTHONFILE) $(CENTROIDS) $(ELEMENTS) > $(KMEANS_INPUT_FILE)
	sed -i 's/INPUT_THREADS/$(NUM_THREADS)/g' kmeans.c

	rm -fr kmeans kmeans.o

	gcc -o kmeans -fopenmp kmeans.c
	for number in 0 1 2; do \
		./kmeans < $(KMEANS_INPUT_FILE) | grep real >> $(KMEANS_OUTPUT_FILE_PAR); \
	done

	rm -fr kmeans kmeans.o

	gcc -o kmeans kmeans.c
	for number in 0 1 2; do \
		{ time ./kmeans < $(KMEANS_INPUT_FILE) >/dev/null; } 2>&1 | grep real >> $(KMEANS_OUTPUT_FILE_SEQ); \
	done

	rm -fr kmeans kmeans.o

	echo "Sequential Num threads: 1, num_elements: $(ELEMENTS), time measured: " >> $(KMEANS_TIME_MEASURE)
	python3 calculate_median.py -f $(KMEANS_OUTPUT_FILE_SEQ) >> $(KMEANS_TIME_MEASURE)
	echo "\n" >> $(KMEANS_TIME_MEASURE)

	echo "Parallel Num threads: $(NUM_THREADS), num_elements: $(ELEMENTS), time measured: " >> $(KMEANS_TIME_MEASURE)
	python3 calculate_median.py -f $(KMEANS_OUTPUT_FILE_PAR) >> $(KMEANS_TIME_MEASURE)
	echo "\n" >> $(KMEANS_TIME_MEASURE)

	rm $(KMEANS_OUTPUT_FILE_SEQ)
	rm $(KMEANS_OUTPUT_FILE_PAR)
	sed -i 's/NUM_THREADS $(NUM_THREADS)/NUM_THREADS INPUT_THREADS/g' kmeans.c

	cat $(KMEANS_TIME_MEASURE)

all : $(clean) $(build_open_mp) $(gen) $(run)
	@echo "Running"
	
# run : $(clean) $(EXE) 
# 	./$(EXE) < $(KMEANS_INPUT_FILE) > $(KMEANS_OUTPUT_FILE)


