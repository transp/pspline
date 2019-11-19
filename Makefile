# Simple makefile to build libpspline.a or libpspline.so

# Include compile options
include Makefile.def

# Define directories and library names
SRCDIR := $(shell pwd)
PSPDIR := $(SRCDIR)/pspline
EZSDIR := $(SRCDIR)/ezspline
CZSDIR := $(SRCDIR)/czspline

PREFIX := $(or $(PREFIX),$(SRCDIR)/build)
LIBDIR := $(PREFIX)/lib
INCDIR := $(PREFIX)/include
TESTDR := $(PREFIX)/test

LIBAR  := libpspline.a
LIBSO  := libpspline.so

EXELUP := lookup_test
EXEPSP := pspltest
EXEQKS := qk_pspline
EXEEZS := ezspline_test
EXECZS := czspline_test

# Define the objects
SRCS :=	$(PSPDIR)/precision_mod.f90 \
	$(filter-out $(PSPDIR)/precision_mod.f90, $(wildcard $(PSPDIR)/*.f90))
SRCS := $(filter-out $(PSPDIR)/$(EXELUP).f90, $(SRCS))
SRCS := $(filter-out $(PSPDIR)/$(EXEPSP).f90, $(SRCS))
SRCS := $(filter-out $(PSPDIR)/$(EXEQKS).f90, $(SRCS))
SRCS += $(EZSDIR)/ezspline_mod.f90 \
	$(filter-out $(EZSDIR)/ezspline_mod.f90, $(wildcard $(EZSDIR)/*.f90))
SRCS := $(filter-out $(EZSDIR)/ezspline_test.f90, $(SRCS))
SRCS += $(wildcard $(CZSDIR)/*.f90)
OBJS :=	$(SRCS:.f90=.o)


.PHONY: clean clobber debug install shared single uninstall

all:	$(LIBAR) $(EXELUP) $(EXEPSP) $(EXEQKS) $(EXEEZS) $(EXECZS)

$(LIBAR): $(OBJS)
	@ar r $@ $(OBJS)

$(EXELUP):
	$(FC) $(FFLAGS) $(PSPDIR)/$@.f90 -o $@ -L$(SRCDIR) -lpspline $(FLIBS)

$(EXEPSP):
	$(FC) $(FFLAGS) $(PSPDIR)/$@.f90 -o $@ -L$(SRCDIR) -lpspline $(FLIBS)

$(EXEQKS):
	$(FC) $(FFLAGS) $(EZSDIR)/$@.f90 -o $@ -L$(SRCDIR) -lpspline $(FLIBS)

$(EXEEZS):
	$(FC) $(FFLAGS) $(EZSDIR)/$@.f90 -o $@ -L$(SRCDIR) -lpspline $(FLIBS)

$(EXECZS):
	$(CXX) $(CXXFLAGS) $(CZSDIR)/$@.cc -o $@ -L$(SRCDIR) -lpspline $(CXXLIBS) $(FLIBS)

install:
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(INCDIR) || mkdir -p $(INCDIR)
	@test -d $(TESTDR) || mkdir -p $(TESTDR)
	@cp *.mod $(INCDIR)
	@test ! -f $(LIBAR) || cp $(LIBAR) $(LIBDIR)
	@test ! -f $(LIBSO) || cp $(LIBSO) $(LIBDIR)
	@cp $(PSPDIR)/pspltest.output $(TESTDR)
	@cp $(EZSDIR)/ezspline_test.ref $(TESTDR)
	@cp $(EZSDIR)/qk_pspline.output $(TESTDR)
	@test ! -f $(EXELUP) || cp $(EXELUP) $(TESTDR)
	@test ! -f $(EXEPSP) || cp $(EXEPSP) $(TESTDR)
	@test ! -f $(EXEQKS) || cp $(EXEQKS) $(TESTDR)
	@test ! -f $(EXEEZS) || cp $(EXEEZS) $(TESTDR)

debug:	FFLAGS += $(DFFLAGS)
debug:	LDFLAGS += $(DFFLAGS)
debug:	$(LIBAR)
debug:	$(EXELUP)
debug:	$(EXEPSP)
debug:  $(EXEQKS)
debug:	$(EXEEZS)
debug:	$(EXECZS)

shared:	FFLAGS += $(SFFLAGS)
shared:	LDFLAGS += $(SFFLAGS)
shared:	$(OBJS)
	$(FC) $(LDFLAGS) -o $(LIBSO) $(OBJS)

single:	FFLAGS += -D_SINGLE_PRECISION
single:	LDFLAGS += -D_SINGLE_PRECISION
single: $(LIBAR)
single: $(EXELUP)
single: $(EXEPSP)
single: $(EXEQKS)
single:	$(EXEEZS)

$(OBJS) :
#
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@ $(FLIBS) $(FINCL)

uninstall:
	@rm -rf $(PREFIX)

clean:
	@rm -f *~ $(PSPDIR)/*~ $(EZSDIR)/*~ $(CZSDIR)/*~
	@rm -f $(OBJS)
	@rm -f *.mod

clobber: clean
	@rm -f $(LIBAR) $(LIBSO)
	@rm -f $(EXELUP) $(EXEPSP) $(EXEQKS) $(EXEEZS) $(EXECZS)
