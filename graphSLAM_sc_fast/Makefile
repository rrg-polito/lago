CXX = g++
CXXFLAGS = -O3 -Wall -I../CSparse/Include
LDFLAGS = -L../CSparse/Lib/ -lcsparse -lm

OBJECTS = main.o
EXECUTABLE = graphSLAM_sc

all:$(EXECUTABLE)

$(EXECUTABLE):$(OBJECTS)
	$(CXX) -o $(EXECUTABLE) $(OBJECTS) $(CXXFLAGS) $(LDFLAGS)   

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)	


	rm *.txt


	


 	
