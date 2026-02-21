CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?=

SRC_DIR := src
BIN_DIR := bin

ATOM_SRC := $(SRC_DIR)/atom.cpp
REALTIME_SRC := $(SRC_DIR)/atom_realtime.cpp
RAYTRACER_SRC := $(SRC_DIR)/atom_raytracer.cpp
WAVE2D_SRC := $(SRC_DIR)/wave_atom_2d.cpp

UNAME_S := $(shell uname -s)

ifeq ($(OS),Windows_NT)
    EXE := .exe
    LDLIBS := -lglfw3 -lglew32 -lopengl32 -lgdi32
else ifeq ($(UNAME_S),Darwin)
    EXE :=
    BREW_PREFIX := $(shell command -v brew >/dev/null 2>&1 && brew --prefix)
    ifneq ($(strip $(BREW_PREFIX)),)
        CXXFLAGS += -I$(BREW_PREFIX)/include
        LDFLAGS += -L$(BREW_PREFIX)/lib
    endif
    LDLIBS := -lglfw -lGLEW -framework OpenGL
else
    EXE :=
    LDLIBS := -lglfw -lGLEW -lGL -lGLU -ldl -lpthread -lm
endif

ATOM_BIN := $(BIN_DIR)/atom$(EXE)
REALTIME_BIN := $(BIN_DIR)/atom_realtime$(EXE)
RAYTRACER_BIN := $(BIN_DIR)/atom_raytracer$(EXE)
WAVE2D_BIN := $(BIN_DIR)/wave_atom_2d$(EXE)

.PHONY: all build atom realtime raytracer wave2d clean

all: build
build: atom realtime raytracer wave2d

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

atom: $(ATOM_BIN)
$(ATOM_BIN): $(ATOM_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

realtime: $(REALTIME_BIN)
$(REALTIME_BIN): $(REALTIME_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

raytracer: $(RAYTRACER_BIN)
$(RAYTRACER_BIN): $(RAYTRACER_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

wave2d: $(WAVE2D_BIN)
$(WAVE2D_BIN): $(WAVE2D_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(ATOM_BIN) $(REALTIME_BIN) $(RAYTRACER_BIN) $(WAVE2D_BIN)
