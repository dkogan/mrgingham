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

$(addsuffix .o,$(basename $(LIB_SOURCES))): CCXXFLAGS += -fvisibility=hidden

CFLAGS    += -std=gnu99
CCXXFLAGS += -Wno-unused-function -Wno-missing-field-initializers -Wno-unused-parameter -Wno-strict-aliasing -Wno-int-to-pointer-cast

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





############# This is a partial copy of mrbuild. Should set this up as a
############# dependency somehow... mrbuild is copyright 2016-2018 California
############# Institute of Technology. Released under the GNU Lesser General
############# Public License version 2.1 or greater.

# This is a common Makefile that can be used as the core buildsystem for
# projects providing a library and some executables using this library. Please
# see README.build.org for the documentation.


# There are two ways to pass variables to make:
#
# make CFLAGS=-foo
#   and
# CFLAGS=-foo make
#
# The former creates a "command line" variable and the latter an
# "environment variable". In order to be able to modify a "command line"
# variable (to add other flags, say), one MUST use 'override'. So one would have to do
#
# override CFLAGS += -bar
#
# without the "override" nothing would happen. I want to avoid this rabbithole
# entirely, so I disallow "command line" variables for things that I modify.
#
# I only do this for xxxFLAGS becuase cross-building and installing in Debian
# uses Make in this way. Hopefully this is safe enough
$(foreach v,$(filter %FLAGS,$(.VARIABLES)),$(if $(filter command line,$(origin $v)), $(error '$v' not allowed as a make parameter. Please do "$v=xxxx make yyyy" instead of "make yyyy $v=xxxx")))



# Make sure I have the variables that must be defined. Libraries need an ABI and
# TAIL version, while other need a plain VERSION
MUST_DEF_VARIABLES        := PROJECT_NAME $(if $(LIB_SOURCES),ABI_VERSION TAIL_VERSION,VERSION)
$(foreach v,$(MUST_DEF_VARIABLES),$(if $($v),,$(error User MUST specify $v)))



# The default VERSION string that appears as a #define to each source file, and
# to any generated documentation (gengetopt and so on). The user can set this to
# whatever they like
VERSION ?= $(ABI_VERSION).$(TAIL_VERSION)

# Default compilers. By default, we use g++ as a linker
CC        ?= gcc
CXX       ?= g++
NVCC      ?= nvcc
CC_LINKER ?= $(CXX)

# used to make gcc output header dependency information. All source
# files generated .d dependency definitions that are included at the
# bottom of this file
CCXXFLAGS += -MMD -MP

# always building with debug information. This is stripped into the
# -dbg/-debuginfo packages by debhelper/rpm later
CCXXFLAGS += -g

# I want the frame pointer. Makes walking the stack WAY easier
CCXXFLAGS += -fno-omit-frame-pointer

# I look through my LIB_SOURCES and BIN_SOURCES. Anything that isn't a wildcard
# (has * or ?) should exist. If it doesn't, the user messed up and I flag it
get_no_wildcards          = $(foreach v,$1,$(if $(findstring ?,$v)$(findstring *,$v),,$v))
complain_if_nonempty      = $(if $(strip $1),$(error $2: $1))
complain_unless_all_exist = $(call complain_if_nonempty,$(call get_no_wildcards,$(filter-out $(wildcard $1),$1)),File not found: )
$(call complain_unless_all_exist,$(LIB_SOURCES) $(BIN_SOURCES))


LIB_SOURCES := $(wildcard $(LIB_SOURCES))
BIN_SOURCES := $(wildcard $(BIN_SOURCES))

LIB_OBJECTS := $(addsuffix .o,$(basename $(LIB_SOURCES)))
BIN_OBJECTS := $(addsuffix .o,$(basename $(BIN_SOURCES)))

SOURCE_DIRS := $(sort ./ $(dir $(LIB_SOURCES) $(BIN_SOURCES)))

# if the PROJECT_NAME is libxxx then LIB_NAME is libxxx
# if the PROJECT_NAME is xxx    then LIB_NAME is libxxx
LIB_NAME           := $(or $(filter lib%,$(PROJECT_NAME)),lib$(PROJECT_NAME))
LIB_TARGET_SO_BARE := $(LIB_NAME).so
LIB_TARGET_SO_ABI  := $(LIB_TARGET_SO_BARE).$(ABI_VERSION)
LIB_TARGET_SO_FULL := $(LIB_TARGET_SO_ABI).$(TAIL_VERSION)
LIB_TARGET_SO_ALL  := $(LIB_TARGET_SO_BARE) $(LIB_TARGET_SO_ABI) $(LIB_TARGET_SO_FULL)

BIN_TARGETS := $(basename $(BIN_SOURCES))


# all objects built for inclusion in shared libraries get -fPIC. We don't build
# static libraries, so this is 100% correct
$(LIB_OBJECTS): CCXXFLAGS += -fPIC

CCXXFLAGS += -DVERSION='"$(VERSION)"'




# These are here to process the options separately for each file being built.
# This allows per-target options to be set
#
# if no explicit optimization flags are given, optimize
define massageopts
$1 $(if $(filter -O%,$1),,-O3)
endef

# If no C++ standard requested, I default to c++0x
define massageopts_cxx
$(call massageopts,$1 $(if $(filter -std=%,$1),,-std=c++0x))
endef

define massageopts_c
$(call massageopts,$1)
endef



# define the compile rules. I need to redefine the rules here because my
# C..FLAGS variables are simple (immediately evaluated), but the user
# could have specified per-target flags that ALWAYS evaluate deferred-ly
#
# I add the warning options AT THE START of the flag list so that the user can
# override these
cc_build_rule = $(strip $(CXX)  $(call massageopts_cxx,-Wall -Wextra $(CXXFLAGS) $(CCXXFLAGS) $(CPPFLAGS))) -c -o $@ $<
c_build_rule  = $(strip $(CC)   $(call massageopts_c,  -Wall -Wextra $(CFLAGS)   $(CCXXFLAGS) $(CPPFLAGS))) -c -o $@ $<
cu_build_rule = $(strip $(NVCC) $(call massageopts_c,  -Wall -Wextra $(CUFLAGS)               $(CPPFLAGS))) -c -o $@ $<
cu_build_rule += --compiler-options "-Wall -Wextra $(CCXXFLAGS)"


%.o:%.C
	$(cc_build_rule)
%.o:%.cc
	$(cc_build_rule)
%.o:%.cpp
	$(cc_build_rule)
%.o:%.c
	$(c_build_rule)
%.o:%.cu
	$(cu_build_rule)
%.o: %.S
	$(CC) $(ASFLAGS) $(CPPFLAGS) -c -o $@ $<

%.h %.c: %.ggo
	gengetopt -C -u -g $(VERSION) -i $< -F $* args -f $(notdir $*) -a $(notdir $*)






# by default I build shared libraries only. We known how to build static
# libraries too, but I don't do it unless asked
all: $(if $(strip $(LIB_SOURCES)),$(LIB_TARGET_SO_ALL)) $(if $(strip $(BIN_SOURCES)),$(BIN_TARGETS))
.PHONY: all
.DEFAULT_GOAL := all

$(LIB_TARGET_SO_FULL): LDFLAGS += -shared -Wl,--default-symver -fPIC -Wl,-soname,$(notdir $(LIB_TARGET_SO_BARE)).$(ABI_VERSION)


$(LIB_TARGET_SO_BARE) $(LIB_TARGET_SO_ABI): $(LIB_TARGET_SO_FULL)
	ln -fs $(notdir $(LIB_TARGET_SO_FULL)) $@

# Here instead of specifying $^, I do just the %.o parts and then the
# others. This is required to make the linker happy to see the dependent
# objects first and the dependency objects last. Same as for BIN_TARGETS
$(LIB_TARGET_SO_FULL): $(LIB_OBJECTS)
	$(CC_LINKER) $(LDFLAGS) $(filter %.o, $^) $(filter-out %.o, $^) $(LDLIBS) -o $@

# I make sure to give the .o to the linker before the .so and everything else.
# The .o may depend on the other stuff. The binaries get an rpath (removed at
# install time)
$(BIN_TARGETS): %: %.o
	$(CC_LINKER) -Wl,-rpath=$(abspath .) $(LDFLAGS) $(filter %.o, $^) $(filter-out %.o, $^) $(LDLIBS) -o $@

# The binaries link with the DSO, if there is one. I need the libxxx.so to build
# the binary, and I need the libxxx.so.abi to run it.
$(BIN_TARGETS): $(if $(strip $(LIB_SOURCES)),$(LIB_TARGET_SO_BARE) $(LIB_TARGET_SO_ABI))




clean:
	rm -rf $(foreach d,$(SOURCE_DIRS),$(addprefix $d,*.a *.o *.so *.so.* *.d)) $(BIN_TARGETS) $(foreach s,.c .h,$(addsuffix $s,$(basename $(shell find . -name '*.ggo')))) $(EXTRA_CLEAN)
distclean: clean

.PHONY: distclean clean



########################### installation business

ifneq (,$(filter install,$(MAKECMDGOALS)))
  ifeq  ($(strip $(DESTDIR)),)
    $(error Tried to make install without having DESTDIR defined \
"make install" is ONLY for package building. \
What are you trying to do?)
  endif
endif


# I process the simple wildcard exceptions on DIST_BIN and DIST_INCLUDE in a
# deferred fashion. The reason is that I wand $(wildcard) to run at install
# time, i.e. after stuff is built, and the files $(wildcard) is looking at
# already exist
DIST_BIN_ORIG     := $(DIST_BIN)
DIST_INCLUDE_ORIG := $(DIST_INCLUDE)

DIST_BIN     = $(filter-out $(wildcard $(DIST_BIN_EXCEPT)),		\
                  $(wildcard $(or $(DIST_BIN_ORIG),    $(BIN_TARGETS))))
DIST_INCLUDE = $(filter-out $(wildcard $(DIST_INCLUDE_EXCEPT)),		\
		  $(wildcard $(DIST_INCLUDE_ORIG)))

ifneq (,$(shell grep -qi debian /etc/os-release 2>/dev/null && echo yep))
  # we're a debian box, use the multiarch dir
  DEB_HOST_MULTIARCH := $(shell dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null)
  USRLIB             := usr/lib/$(DEB_HOST_MULTIARCH)
else
  # we're something else. If /usr/lib64 exists, use that. Otherwise /usr/lib
  USRLIB := $(if $(wildcard /usr/lib64),usr/lib64,usr/lib)
endif

# Generates the install rules. Arguments:
#   1. variable containing the being installed
#   2. target path they're being installed to
#   3. post-install commands
define install_rule
$(if $(strip $($1)),
	mkdir -p $2 &&									\
	cp -r $($1) $2 &&								\
	$(if $($(1)_EXCEPT_FINDSPEC),find $2 $($(1)_EXCEPT_FINDSPEC) -delete &&)	\
	$(or $3,true) )
endef

ifneq ($(strip $(LIB_SOURCES)),)
install: $(LIB_TARGET_SO_ALL)
endif

install: $(BIN_TARGETS) $(DIST_DOC) $(DIST_MAN) $(DIST_DATA)

# using 'cp -P' instead of 'install' because the latter follows links unconditionally
ifneq ($(strip $(LIB_SOURCES)),)
	mkdir -p $(DESTDIR)/$(USRLIB)
	cp -P $(LIB_TARGET_SO_FULL)  $(DESTDIR)/$(USRLIB)
	ln -fs $(notdir $(LIB_TARGET_SO_FULL)) $(DESTDIR)/$(USRLIB)/$(notdir $(LIB_TARGET_SO_ABI))
	ln -fs $(notdir $(LIB_TARGET_SO_FULL)) $(DESTDIR)/$(USRLIB)/$(notdir $(LIB_TARGET_SO_BARE))
endif
	$(call install_rule,DIST_BIN,         $(DESTDIR)/usr/bin,)
	$(call install_rule,DIST_INCLUDE,     $(DESTDIR)/usr/include/$(PROJECT_NAME),)
	$(call install_rule,DIST_DOC,         $(DESTDIR)/usr/share/doc/$(PROJECT_NAME),)
	$(call install_rule,DIST_MAN,         $(DESTDIR)/usr/share/man,)
	$(call install_rule,DIST_DATA,        $(DESTDIR)/usr/share/$(PROJECT_NAME),)
	$(call install_rule,DIST_PERL_MODULES,$(DESTDIR)/usr/share/perl5,)
	$(call install_rule,DIST_PY2_MODULES, $(DESTDIR)/usr/lib/python2.7/site-packages,)
	$(call install_rule,DIST_PY3_MODULES, $(DESTDIR)/usr/lib/python3.4/site-packages,)

        # In filenames I rename __colon__ -> :
        # This is required because Make can't deal with : in rules
	for fil in `find $(DESTDIR) -name '*__colon__*'`; do mv $$fil `echo $$fil | sed s/__colon__/:/g`; done

        # Remove rpaths from everything. /usr/bin is allowed to fail because
        # some of those executables aren't ELFs. On the other hand, any .so we
        # find IS en ELF. Those could live in a number of places, since they
        # could be extension modules for the various languages, and I thus look
        # for those EVERYWHERE
ifneq ($(strip $(DIST_BIN)),)
	chrpath -d $(DESTDIR)/usr/bin/* 2>/dev/null || true
endif
	find $(DESTDIR) -name '*.so' | xargs chrpath -d

        # Any perl programs need their binary path stuff stripped out. This
        # exists to let these run in-tree, but needs to be removed at
        # install-time (similar to an RPATH)
ifneq ($(strip $(DIST_BIN)),)
	for fil in `find $(DESTDIR)/usr/bin -type f`; do head -n1 $$fil | grep -q '^#!.*/perl$$' && perl -n -i -e 'print unless /^\s* use \s+ lib \b/x' $$fil || true; done
endif
	test -e $(DESTDIR)/usr/share/perl5 && find $(DESTDIR)/usr/share/perl5 -type f | xargs perl -n -i -e 'print unless /^\s* use \s+ lib \b/x' || true

ifneq ($(strip $(DIST_BIN)),)
        # Everything executable needs specific permission bits
	chmod 0755 $(DESTDIR)/usr/bin/*
endif

.PHONY: install



# I want to keep all the intermediate files always
.SECONDARY:

# the header dependencies
-include $(addsuffix *.d,$(SOURCE_DIRS))

