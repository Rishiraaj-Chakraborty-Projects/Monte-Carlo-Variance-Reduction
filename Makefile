CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS  ?= 
BIN       = mc

all: $(BIN)

$(BIN): src/main.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN) src/main.cpp $(LDFLAGS)

test: $(BIN)
	./$(BIN) --test

run: $(BIN)
	./$(BIN)

clean:
	rm -f $(BIN) $(BIN).exe

.PHONY: all test run clean
