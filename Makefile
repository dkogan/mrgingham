PROJECT_NAME := mrgingham
ABI_VERSION  := 0
TAIL_VERSION := 1

BIN_SOURCES := mrgingham_test_from_image mrgingham_test_from_points
BIN_SOURCES += test_dump_blobs
LIB_SOURCES := find_grid.cc find_blobs.cc mrgingham.cc

CXXFLAGS_CV := $(shell pkg-config --cflags opencv)
LDLIBS_CV   := $(shell pkg-config --libs   opencv)
CCXXFLAGS += $(CXXFLAGS_CV)
LDLIBS    += $(LDLIBS_CV)

CCXXFLAGS += -Wno-unused-function -Wno-missing-field-initializers

include /usr/include/mrbuild/Makefile.common
