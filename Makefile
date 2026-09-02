CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude
LDFLAGS  ?= -static

LIB_SRC   = src/simulator.cpp
DEMO_BIN  = mc_demo
TEST_BIN  = mc_test
BENCH_BIN = mc_bench

all: $(DEMO_BIN) $(TEST_BIN) $(BENCH_BIN)

$(DEMO_BIN): apps/main.cpp $(LIB_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(TEST_BIN): tests/test_monte_carlo.cpp $(LIB_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BENCH_BIN): benchmarks/bench_monte_carlo.cpp $(LIB_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

bench: $(BENCH_BIN)
	./$(BENCH_BIN)

run: $(DEMO_BIN)
	./$(DEMO_BIN)

clean:
	rm -f $(DEMO_BIN) $(DEMO_BIN).exe $(TEST_BIN) $(TEST_BIN).exe $(BENCH_BIN) $(BENCH_BIN).exe

.PHONY: all test bench run clean
