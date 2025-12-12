# ========================================
#   Compiler & Flags
# ========================================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -march=native -fPIC \
            -fno-math-errno -ffast-math \
            -Wall -Wextra -Iinclude

# ========================================
#   GoogleTest
# ========================================
GTEST_ROOT := thirdparty/googletest
GTEST_INC  := -I$(GTEST_ROOT)/googletest/include -I$(GTEST_ROOT)/googletest
CXXFLAGS  += $(GTEST_INC)

$(GTEST_ROOT):
	@echo "GoogleTest not found → downloading..."
	@mkdir -p thirdparty
	@git clone https://github.com/google/googletest.git $(GTEST_ROOT)

GTEST_ALL_SRC := $(GTEST_ROOT)/googletest/src/gtest-all.cc
GTEST_OBJ     := build/gtest-all.o

$(GTEST_OBJ): $(GTEST_ALL_SRC) | $(GTEST_ROOT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ========================================
#   Python / Pybind
# ========================================
PYTHON := python3
PYTHON_INCLUDES := $(shell $(PYTHON)-config --includes)
PYTHON_LDFLAGS  := $(shell $(PYTHON)-config --ldflags | sed 's/-lpython[0-9.]*//')
PYEXT := $(shell $(PYTHON)-config --extension-suffix)

PY_SRC_DIR := src/python
PY_MODULE  := ftst_dense$(PYEXT)
PY_MODULE_PATH := $(PY_SRC_DIR)/$(PY_MODULE)

PYBIND_INCLUDES := -I/usr/include/pybind11

# ========================================
#   Source / Build
# ========================================
SRC_DIR := src/cpp

SRC := $(SRC_DIR)/main.cpp $(SRC_DIR)/solver.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := solver_app

TEST_SRC := $(SRC_DIR)/solver_tests.cpp
TEST_OBJ := $(TEST_SRC:.cpp=.o)
TEST_TARGET := solver_tests

PYBIND_SRC := $(SRC_DIR)/pybind_module.cpp
PYBIND_OBJ := $(PYBIND_SRC:.cpp=.o)

# ========================================
#   Rules
# ========================================
all: $(TARGET) $(TEST_TARGET) $(PY_MODULE_PATH)

# Main solver
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ -fopenmp -lpthread -ldl -lm

# GoogleTest binary
$(TEST_TARGET): $(TEST_OBJ) $(GTEST_OBJ) $(SRC_DIR)/solver.o
	$(CXX) $(TEST_OBJ) $(SRC_DIR)/solver.o $(GTEST_OBJ) \
	    -o $@ -fopenmp -lpthread -ldl -lm

# Pybind object
$(PYBIND_OBJ): $(PYBIND_SRC)
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) $(PYBIND_INCLUDES) -c $< -o $@

# Python extension (.so) — top-level module
$(PY_MODULE_PATH): $(PYBIND_OBJ) $(SRC_DIR)/solver.o
	@mkdir -p $(PY_SRC_DIR)
	$(CXX) -shared $(PYBIND_OBJ) $(SRC_DIR)/solver.o \
	    -o $@ $(PYTHON_LDFLAGS) -fopenmp
	@echo "Built Python extension → $(PY_MODULE_PATH)"

# Generic object rule
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ========================================
#   Utilities
# ========================================
clean:
	rm -f $(OBJ) $(TARGET) \
	      $(TEST_OBJ) $(TEST_TARGET) \
	      $(GTEST_OBJ) \
	      $(PYBIND_OBJ) $(PY_MODULE_PATH)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

pytest: $(PY_MODULE_PATH)
	PYTHONPATH="$(CURDIR)/src/python" pytest -q tests

run: $(TARGET)
	./$(TARGET)

print:
	@echo "Python module path = $(PY_MODULE_PATH)"

.PHONY: all clean run print test pytest
