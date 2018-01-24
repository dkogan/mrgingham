PROJECT_NAME := mrgingham
ABI_VERSION  := 1
TAIL_VERSION := 1


DIST_BIN := mrgingham_test_from_image mrgingham_test_from_points
BIN_SOURCES := $(DIST_BIN:%=%.cc) test_dump_blobs.cc test_dump_chessboard_corners.cc
LIB_SOURCES := find_grid.cc find_blobs.cc find_chessboard_corners.cc mrgingham.cc ChESS.c

CXXFLAGS_CV := $(shell pkg-config --cflags opencv)
LDLIBS_CV   := $(shell pkg-config --libs   opencv)
CCXXFLAGS += $(CXXFLAGS_CV)
LDLIBS    += $(LDLIBS_CV)

$(addsuffix .o,$(basename $(LIB_SOURCES))): CCXXFLAGS += -fvisibility=hidden

CFLAGS    += -std=gnu99
CCXXFLAGS += -Wno-unused-function -Wno-missing-field-initializers -Wno-unused-parameter

DIST_INCLUDE := mrgingham.hh point.hh

include /usr/include/mrbuild/Makefile.common
