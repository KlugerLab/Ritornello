CXXFLAGS = -Wall -fmessage-length=0
INCLUDEPATH = -I include/ -I /usr/include/samtools/ -I ~/fftw3Installation/include/ -I ~/samtools-0.1.19/
LIBS = -L ~/samtools-0.1.19/ -L ~/fftw3Installation/lib/ -lbam -lz -lm #-lblas

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix bin/,$(notdir $(CPP_FILES:.cpp=.o)))
TARGET_FILE =	bin/Ritornello

#builds target
all: CXXFLAGS += -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -fopenmp
all: LIBS += -l:libfftw3.a -l:libfftw3_omp.a
all:	$(TARGET_FILE)

debug: CXXFLAGS += -O0 -DDEBUG -g -fopenmp
debug: LIBS += -l:libfftw3.a -l:libfftw3_omp.a
#debug: LIBS += -lfftw3
debug: $(TARGET_FILE)

#cleans all objects
clean:
	rm -f $(OBJ_FILES) $(TARGET_FILE)

#links all objects into the library
$(TARGET_FILE):	$(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $(TARGET_FILE) $(OBJ_FILES) $(LIBS)
	
#Builds all CPP files into objects	
bin/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEPATH)		
