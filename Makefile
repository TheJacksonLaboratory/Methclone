CXX = g++
CXXOPT = -O3
CPPFLAGS += -Igzstream
CPPFLAGS += -IBamTools/include
LDFLAGS += -Lgzstream
LDFLAGS += -LBamTools/lib 
LDFLAGS += -lbamtools -lgzstream -lz
PROGRAME = methClone



all: $(PROGRAME)


$(PROGRAME): OptionParser.cpp methClone.cpp
	$(CXX) $(CXXOPT) $^ $(CPPFLAGS) $(LDFLAGS) -o $@

clean:
	rm  methClone
