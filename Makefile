# ========================================
#   Compiler & Flags
# ========================================
export PYTHONPATH := $(CURDIR)/src/python:$(PYTHONPATH)
CXX := g++
CXXFLAGS = -O3 -std=c++17 -march=native -fPIC \
           -fno-math-errno -ffast-math \
           -Wall -Wextra

# ========================================
#   GTest Include
# ========================================
GTEST_ROOT := thirdparty/googletest
GTEST_FLAGS := -I$(GTEST_ROOT)/googletest/include -I$(GTEST_ROOT)/googletest
CXXFLAGS += $(GTEST_FLAGS)
CXXFLAGS += -Iinclude
# ========================================
#   GTest Include
# ========================================
GTEST_ROOT := thirdparty/googletest
GTEST_FLAGS := -I$(GTEST_ROOT)/googletest/include -I$(GTEST_ROOT)/googletest
CXXFLAGS += $(GTEST_FLAGS)
CXXFLAGS += -Iinclude


$(GTEST_ROOT):
	@echo "GoogleTest not found → downloading..."
	@mkdir -p thirdparty
	@git clone https://github.com/google/googletest.git $(GTEST_ROOT)
	@echo "GoogleTest downloaded."

GTEST_ALL_SRC := $(GTEST_ROOT)/googletest/src/gtest-all.cc
GTEST_OBJ := build/gtest-all.o

$(GTEST_OBJ): $(GTEST_ALL_SRC) | $(GTEST_ROOT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $(GTEST_ALL_SRC) -o $@

$(TEST_TARGET): $(TEST_OBJ) $(GTEST_OBJ) $(SRC_DIR)/solver.o | $(GTEST_ROOT)
	$(CXX) $(TEST_OBJ) $(SRC_DIR)/solver.o $(GTEST_OBJ) \
	    -o $@ -fopenmp -lpthread -ldl -lm


# ========================================
#   Python & Pybind
# ========================================
PYTHON       := python3
PYTHON_INCLUDES := $(shell $(PYTHON)-config --includes)
PYTHON_LDFLAGS  := $(shell $(PYTHON)-config --ldflags | sed 's/-lpython[0-9.]*//')
PYEXT := $(shell $(PYTHON)-config --extension-suffix)

PY_SRC_DIR := src/python
PY_PKG_DIR := $(PY_SRC_DIR)/FTST
PY_MODULE  := ftst_dense$(PYEXT)
PY_MODULE_PATH := $(PY_PKG_DIR)/$(PY_MODULE)

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
GTEST_OBJ := build/gtest-all.o
TEST_TARGET := solver_tests

PYBIND_SRC := $(SRC_DIR)/pybind_module.cpp
PYBIND_OBJ := $(PYBIND_SRC:.cpp=.o)

# ========================================
#   Rules
# ========================================
all: $(TARGET) $(TEST_TARGET) $(PY_MODULE_PATH)

# C++ Main solver
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ -fopenmp -lpthread -ldl -lm

# GTest build
$(TEST_TARGET): $(TEST_OBJ) $(GTEST_OBJ) $(SRC_DIR)/solver.o
	$(CXX) $(TEST_OBJ) $(SRC_DIR)/solver.o $(GTEST_OBJ) \
	    -o $@ -fopenmp -lpthread -ldl -lm

# Pybind object
$(PYBIND_OBJ): $(PYBIND_SRC)
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) $(PYBIND_INCLUDES) -c $< -o $@

# Python module output (.so)
$(PY_MODULE_PATH): $(PYBIND_OBJ) $(SRC_DIR)/solver.o
	@mkdir -p $(PY_PKG_DIR)
	$(CXX) -shared $(PYBIND_OBJ) $(SRC_DIR)/solver.o \
	    -o $(PY_MODULE_PATH) \
	    $(PYTHON_LDFLAGS) -fopenmp
	@echo "Built Python extension → $(PY_MODULE_PATH)"

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# GTest all-in-one
$(GTEST_OBJ): $(GTEST_ROOT)/googletest/src/gtest-all.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ========================================
#   Python Packaging (virtualenv)
# ========================================
PYPROJECT := pyproject.toml
VENV_DIR := .venv
VENV_PIP := $(VENV_DIR)/bin/pip

# Auto-generate pyproject.toml if missing
$(PYPROJECT):
	@echo "[build-system]" > $(PYPROJECT)
	@echo "requires = [\"setuptools\", \"wheel\"]" >> $(PYPROJECT)
	@echo "build-backend = \"setuptools.build_meta\"" >> $(PYPROJECT)
	@echo "" >> $(PYPROJECT)
	@echo "[project]" >> $(PYPROJECT)
	@echo "name = \"FTST\"" >> $(PYPROJECT)
	@echo "version = \"0.1.0\"" >> $(PYPROJECT)
	@echo "packages = [\"FTST\"]" >> $(PYPROJECT)
	@echo "package-dir = {\"\"=\"src/python\"}" >> $(PYPROJECT)
	@echo "" >> $(PYPROJECT)
	@echo "[tool.setuptools.package-data]" >> $(PYPROJECT)
	@echo "FTST = [\"*.so\"]" >> $(PYPROJECT)
	@echo "Generated pyproject.toml."

$(VENV_DIR):
	$(PYTHON) -m venv $(VENV_DIR)
	$(VENV_PIP) install --upgrade pip setuptools wheel

install: $(PYPROJECT) $(PY_MODULE_PATH) $(VENV_DIR)
	$(VENV_PIP) install .
	@echo ""
	@echo "========================================"
	@echo " FTST installed inside virtualenv: .venv"
	@echo " Activate with:"
	@echo "   source .venv/bin/activate"
	@echo "========================================"

uninstall:
	@$(VENV_PIP) uninstall -y FTST

wheel: $(PYPROJECT) $(PY_MODULE_PATH) $(VENV_DIR)
	$(VENV_PIP) wheel . -w dist
	@echo "Wheel built → dist/"

# ========================================
#   Utilities
# ========================================
clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_TARGET) \
	      $(GTEST_OBJ) $(PYBIND_OBJ) $(PY_MODULE_PATH)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

pytest: $(PY_MODULE_PATH)
	PYTHONPATH="$(CURDIR)/$(PY_SRC_DIR)" pytest -q tests

run: $(TARGET)
	./$(TARGET)

print:
	@echo "Python module path = $(PY_MODULE_PATH)"

.PHONY: all clean run print test pytest install uninstall wheel
