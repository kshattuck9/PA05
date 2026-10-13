CXX = $(shell brew --prefix llvm)/bin/clang++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp
CPPFLAGS = -Iinclude -Iexamples -I$(OMP_PREFIX)/include
IMPL = src/sparse.cpp
OMP_PREFIX = $(shell brew --prefix libomp)
.PHONY: all test submission clean
CXX = $(shell brew --prefix llvm)/bin/clang++
OMP_PREFIX = $(shell brew --prefix libomp)
LDFLAGS = -L$(OMP_PREFIX)/lib -Wl,-rpath,$(OMP_PREFIX)/lib -lomp
all: build/demo build/bench build/check build/refinement build/compare
build:
	mkdir -p build
build/sparse.o: $(IMPL) include/sparse.h | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(IMPL) -o $@
build/libsparse.a: build/sparse.o
	ar rcs $@ $^
build/demo: examples/demo.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/bench: examples/bench.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/check: tests/check.cpp build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a $(LDFLAGS) -o $@
test: build/check
	./build/check
submission: | build
	python3 -c 'import zipfile; z=zipfile.ZipFile("build/submission.zip","w",zipfile.ZIP_DEFLATED); z.write("src/sparse.cpp","sparse.cpp"); z.close()'
clean:
	rm -rf build
build/refinement: examples/refinement.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/compare: examples/compare.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/parallel: examples/parallel.cpp build/libsparse.a | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) examples/parallel.cpp build/libsparse.a $(LDFLAGS) -o $@