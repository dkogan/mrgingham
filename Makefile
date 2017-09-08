PROJECT_NAME := mrgrid
ABI_VERSION  := 0
TAIL_VERSION := 0

BIN_SOURCES := test_from_image.cc test_from_points.cc
LIB_SOURCES := find_grid.cc find_blobs.cc

CXXFLAGS_CV := $(shell pkg-config --cflags opencv)
LDLIBS_CV   := $(shell pkg-config --libs   opencv)
CCXXFLAGS += $(CXXFLAGS_CV)
LDLIBS    += $(LDLIBS_CV)

CCXXFLAGS += -Wno-unused-function -Wno-missing-field-initializers

include /usr/include/mrbuild/Makefile.common
