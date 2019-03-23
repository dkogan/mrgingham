include Makefile.common.header


PROJECT_NAME := mrgingham
ABI_VERSION  := 1
TAIL_VERSION := 1


DIST_BIN := mrgingham mrgingham-observe-pixel-uncertainty
DIST_MAN := $(addsuffix .1,$(DIST_BIN))

# I want the tool I ship to be called "mrgingham", but I already have
# mrgingham.cc: it's a part of the LIBRARY
mrgingham: mrgingham-from-image
	cp $< $@
EXTRA_CLEAN += mrgingham
all: mrgingham

BIN_SOURCES := mrgingham-from-image.cc
BIN_SOURCES += test-dump-chessboard-corners.cc test-dump-blobs.cc test-find-grid-from-points.cc

LIB_SOURCES := find_grid.cc find_blobs.cc find_chessboard_corners.cc mrgingham.cc ChESS.c

CXXFLAGS_CV := $(shell pkg-config --cflags opencv)
LDLIBS_CV   := $(shell pkg-config --libs   opencv)
CCXXFLAGS += $(CXXFLAGS_CV)
LDLIBS    += $(LDLIBS_CV) -lpthread

CCXXFLAGS += -fvisibility=hidden

CFLAGS    += -std=gnu99
CCXXFLAGS += -Wno-unused-function -Wno-missing-field-initializers -Wno-unused-parameter -Wno-strict-aliasing -Wno-int-to-pointer-cast -Wno-unused-variable

# opencv is a ridiculous library that keeps breaking their APIs for some reason.
# On opencv2 I do this:
AM_OLD_RH7 := $(shell grep "release 7" /etc/redhat-release 2>/dev/null)
CCXXFLAGS  += -DOLD_OPENCV=$(if $(AM_OLD_RH7),1,0)

# On opencv4 I do this:
ifneq ($(wildcard /usr/include/opencv4),)
CCXXFLAGS += -D CV_LOAD_IMAGE_GRAYSCALE=cv::IMREAD_GRAYSCALE
endif


DIST_INCLUDE := mrgingham.hh point.hh

# I construct the README.org from the template. The only thing I do is to insert
# the manpages. Note that this is more complicated than it looks:
#
# 1. The documentation lives in a POD
# 2. This documentation is stripped out here with pod2text, and included in the
#    README. This README is an org-mode file, and the README.template.org
#    container included the manpage text inside a #+BEGIN_EXAMPLE/#+END_EXAMPLE.
#    So the manpages are treated as a verbatim, unformatted text blob
# 3. Further down, the same POD is converted to a manpage via pod2man
define MAKE_README =
BEGIN									\
{									\
  for $$a (@ARGV)							\
  {									\
    $$base = $$a =~ s/\.pod$$//r;                                       \
    $$c{$$base} = `pod2text $$a | mawk "/REPOSITORY/{exit} {print}"`;	\
  }									\
}									\
									\
while(<STDIN>)								\
{									\
  print s/xxx-manpage-(.*?)-xxx/$$c{$$1}/gr;				\
}
endef

README.org: README.template.org $(DIST_BIN:%=%.pod)
	< $(filter README%,$^) perl -e '$(MAKE_README)' $(filter-out README%,$^) > $@
all: README.org

%.1: %.pod
	pod2man --center="mrgingham: chessboard corner finder" --name=MRGINGHAM --release="mrgingham $(VERSION)" --section=1 $^ $@
mrgingham-observe-pixel-uncertainty.pod: %.pod: %
	./make-pod.pl $< > $@
	cat footer.pod >> $@
EXTRA_CLEAN += *.1 mrgingham-observe-pixel-uncertainty.pod README.org

########## python stuff

# In the python api I have to cast a PyCFunctionWithKeywords to a PyCFunction,
# and the compiler complains. But that's how Python does it! So I tell the
# compiler to chill
mrgingham_pywrap.o: CFLAGS += -Wno-cast-function-type

mrgingham_pywrap.o mrgingham_pywrap_cplusplus_bridge.o: CCXXFLAGS += $(PY_MRBUILD_CFLAGS)
mrgingham_pywrap.o: $(addsuffix .h,$(wildcard *.docstring))

# The python library is called "mrgingham.so". This is confusing, but is the
# best I could do. Because I want to be able to "import mrgingham"; and the
# normal way of creating a "mrgingham" subdirectory for all the python stuff
# doesn't work here: I already have a directory entry called "mrgingham"; it's
# the main commandline tool.
mrgingham$(PY_EXT_SUFFIX): mrgingham_pywrap.o mrgingham_pywrap_cplusplus_bridge.o libmrgingham.so
	$(PY_MRBUILD_LINKER) $(PY_MRBUILD_LDFLAGS) $^ -o $@

DIST_PY2_MODULES := mrgingham$(PY_EXT_SUFFIX)

all: mrgingham$(PY_EXT_SUFFIX)

include Makefile.common.footer
