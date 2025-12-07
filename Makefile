# ============================
#   Compiler & Flags
# ============================
CXX := g++
CXXFLAGS := -O3 -std=gnu++17 -Wall -Wextra -Wpedantic -fPIC -fopenmp

GTEST_ROOT := ../FTST/thirdparty/googletest
GTEST_FLAGS := -I$(GTEST_ROOT)/googletest/include -I$(GTEST_ROOT)/googletest
CXXFLAGS += $(GTEST_FLAGS)

PYTHON_INCLUDES := $(shell python3-config --includes)
PYTHON_LDFLAGS := $(shell python3-config --ldflags)
PYEXT := $(shell python3-config --extension-suffix)
PYBIND_INCLUDES := -I/usr/include/pybind11
PY_MODULE := ftst_dense$(PYEXT)

# ============================
#   Source / Build
# ============================
SRC := main.cpp solver.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := solver_app

TEST_SRC := tests/solver_tests.cpp
TEST_OBJ := $(TEST_SRC:.cpp=.o)
GTEST_OBJ := build/gtest-all.o
TEST_TARGET := solver_tests

PYBIND_SRC := pybind_module.cpp
PYBIND_OBJ := $(PYBIND_SRC:.cpp=.o)

# ============================
#   Build Rules
# ============================
all: $(TARGET) $(TEST_TARGET) $(PY_MODULE)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ -fopenmp -lpthread -ldl -lm

$(TEST_TARGET): $(TEST_OBJ) $(GTEST_OBJ) solver.o
	$(CXX) $(TEST_OBJ) solver.o $(GTEST_OBJ) -o $@ -fopenmp -lpthread -ldl -lm

$(PYBIND_OBJ): $(PYBIND_SRC)
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) $(PYBIND_INCLUDES) -c $< -o $@

$(PY_MODULE): $(PYBIND_OBJ) solver.o
	$(CXX) -shared $(PYBIND_OBJ) solver.o -o $@ $(PYTHON_LDFLAGS) -fopenmp

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(GTEST_OBJ): $(GTEST_ROOT)/googletest/src/gtest-all.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ============================
#   Utilities
# ============================
clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_TARGET) $(GTEST_OBJ) $(PYBIND_OBJ) $(PY_MODULE)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

pytest: $(PY_MODULE)
	PYTHONPATH="$(CURDIR)" pytest -q tests_py

run: $(TARGET)
	./$(TARGET)

print:
	@echo "CXX = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "Target = $(TARGET)"

.PHONY: all clean run print test pytest
