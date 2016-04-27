.PHONY : env clean

CC = g++
CCF = gfortran
AR = ar
C_FLAGS = -O2  -Wall -fPIC 

ROOTSYS = $(shell root-config --prefix)
INCLUDE = -I./ -I$(ROOTSYS)/include

CERN_LIBS = -L/usr/lib64/cernlib/2006/lib \
-lgeant321 -lpawlib -lgraflib -lgrafX11 \
-lpacklib -lmathlib -lkernlib 

ROOT_LIBS = $(shell root-config --libs) 
ROOT_LIBS += -lEG -lEGPythia6

MISC_LIBS = ./libbggenXsec.a

LIBS = $(MISC_LIBS) $(CERN_LIBS) $(ROOT_LIBS)

SOURCES = $(shell ls *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = vsGen

SOURCESF = $(shell ls *.F)
OBJECTSF = $(SOURCESF:.F=.o)

all: libbggenXsec.a  $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) libbggenXsec.a
	$(CC) $(C_FLAGS) $^ -o $@ $(LIBS)

libbggenXsec.a: $(OBJECTSF)
	$(AR) r libbggenXsec.a $(OBJECTSF)

#libbggenXsec.a:  $(OBJECTSF)

%.o: %.F
	$(CCF) -c $(C_FLAGS) $(INCLUDE) $<

%.o: %.cpp
	$(CC) -c $(C_FLAGS) $(INCLUDE) $<

clean:
	rm -rf *.o *.a $(EXECUTABLE)

env:
	@echo CC = $(CC)
	@echo CCF = $(CCF)
	@echo C_FLAGS = $(C_FLAGS)
	@echo MAIN_FILE = $(MAIN_FILE)
	@echo EXECUTABLE = $(EXECUTABLE)
	@echo SOURCES = $(SOURCES)
	@echo SOURCESF = $(SOURCESF)
	@echo OBJECTS = $(OBJECTS)
	@echo OBJECTSF = $(OBJECTSF)
	@echo INCLUDE = $(INCLUDE)
	@echo LIBS = $(LIBS)
