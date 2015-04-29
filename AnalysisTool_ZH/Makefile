CC = g++ 
CFLAGS = -O2 -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare
INCLUDES = -I $(shell root-config --incdir) -I$(ATINCLUDE) 
ROOTLIBS = $(shell root-config --libs) -lMinuit -lMathCore -lrt -L$(ATLIB) -lAnalysisTool 
OBJ = main.o 
 
myanalysis: $(OBJ) 
	$(CC) $(CFLAGS) $(ROOTLIBS) -o myanalysis $(OBJ) 

%.o: %.cc 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 
