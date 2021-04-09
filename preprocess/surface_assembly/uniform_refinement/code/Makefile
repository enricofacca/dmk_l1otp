include Rules_make/Rules.make

OBJDIR        = $(CURDIR)/$(ODIR)
MODDIR        = $(CURDIR)/$(MDIR)
LIBDIR        = $(CURDIR)/$(LDIR)
DOXDIR        = $(CURDIR)/$(DDIR)
DOCDIR        = $(CURDIR)/doc
RUNDIR        = $(CURDIR)/run1

INCLUDE       += -I$(MODDIR)

F90FLAGS      = $(OFLAGS) $(FLAGSF90) $(MODEL)  $(INCLUDE) 

FFLAGS        = $(OFLAGS) $(MODEL)  $(INCLUDE) 

LIBSRCS       = 

DEST	      = .

EXTHDRS	      =

HDRS	      = 

LIBS	      = -lblas -llapack

LDFLAGS	      = $(MODEL)  $(INCLUDE) -L. -L/usr/local/lib -L$(LIBDIR)

LINKER	      = $(FC)

MAKEFILE      = Makefile

PRINT	      = pr

CAT	      = cat

PROGRAM	      = uniform_refinement.out

SRCS	      = 00modGlobals.f90 10modMatrix.f90\
		22modSparseMatrix.f90 50modGeometry.f90 main.f90
OBJS	      = 00modGlobals.o 10modMatrix.o 22modSparseMatrix.o\
		50modGeometry.o main.o
SRCF90       := $(filter %.f90,$(SRCS))
SRCF         := $(filter %.f,$(SRCS))
OBJECTS90 = $(SRCF90:%.f90=$(OBJDIR)/%.o) 
OBJECTSF  = $(SRCF:%.f=$(OBJDIR)/%.o) 
OBJECTS   = $(OBJECTS90) $(OBJECTSF) 

PRJS= $(SRCS:jo=.prj)

MODULES        = $(addprefix $(MODDIR)/,$(MODS))

.SUFFIXES: .prj .f90

print-%  : 
	@echo $* = $($*)

.f.prj:
	ftnchek -project -declare -noverbose $<

.f90.o:
	$(FC) $(F90FLAGS) $(INCLUDE) -c  $< 

.f.o:
	$(FC) $(FFLAGS) $(INCLUDE) -c  $< 

#all::
#		set -e; for i in $(SUBDIRS); do cd $$i; $(MAKE); cd ..; done

all::		$(PROGRAM) 

$(PROGRAM):     $(LIBS) $(OBJECTS)
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJECTS) $(LIBS)

$(LIBS):
		set -e; for i in $(LIBSRCS); do cd $$i; $(MAKE); cd ..; done


#$(OBJECTSF90): $(OBJDIR)/%.o: %.f90 
$(OBJDIR)/%.o: %.f90 
		$(FC) $(CPPFLAGS) $(F90FLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.f 
		$(FC) $(CPPFLAGS) $(FFLAGS) -o $@ -c $<

.PHONY:         doc

doc: 
		@-rm -f $(DOXDIR)/README.md
		@cat $(DOXDIR)/README.tpl $(DOXDIR)/*.txt > $(DOXDIR)/README.md
		doxygen $(DOXDIR)/DoxygenConfigFortran 

dirs: 
		@-mkdir -p $(DOCDIR) $(OBJDIR) $(MODDIR) $(LIBDIR)

clean-emacs:
		@-rm -f $(CURDIR)/*.*~ 
		@-rm -f $(CURDIR)/#*# 

check: $(PRJS)
	ftnchek -noverbose -declare $(PRJS) -project -noextern -library > $(PROGRAM).ftn

profile:;       @make "F90FLAGS=$(PFLAGS) $(MODEL) " "FFLAGS=$(PFLAGS) $(MODEL) " "CFLAGS=$(PFLAGS) $(MODEL)" "LDFLAGS=$(PFLAGS) $(LDFLAGS)" $(PROGRAM)

debug:;         @make "F90FLAGS=$(DFLAGS) $(MODEL) $(INCLUDE)" "FFLAGS=$(DFLAGS) $(MODEL) $(INCLUDE)" "LDFLAGS=$(DFLAGS) $(LDFLAGS)" $(PROGRAM)

openmp:;         @make "F90FLAGS=$(OFLAGS) $(OPENMP) $(MODEL) $(INCLUDE)" "FFLAGS=$(OFLAGS) $(OPENMP) $(MODEL) $(INCLUDE)" "LDFLAGS=$(LDFLAGS) $(OPENMP)" $(PROGRAM)

clean:;		@rm -f $(OBJECTS) $(MODULES) $(PROGRAM).cat $(PROGRAM).ftn
		@set -e; for i in $(LIBSRCS); do cd $$i; $(MAKE) clean; cd ..; done

clobber:;	@rm -f $(OBJECTS) $(MODULES) $(PROGRAM).cat $(PROGRAM).ftn $(PROGRAM)
		@-rm -rf $(DOCDIR)/latex/* $(DOCDIR)/html/*
		@-rm -f $(CURDIR)/*.*~ 
		@-rm -f $(RUNDIR)/output/* $(RUNDIR)/output/vtkout/*
		@set -e; for i in $(LIBSRCS); do cd $$i; $(MAKE) clobber; cd ..; done

#depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)
#depend:;	./scripts/mkdeps >> $(MAKEFILE)
#depend:;	./scripts/mkdeps > Rules_make/dependencies

.PHONY:		mods

mods:;		./scripts/find_mod_names >> $(MAKEFILE)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

print:;		$(PRINT) $(HDRS) $(SRCS)

cat:;		$(CAT) $(HDRS) $(SRCS) > $(PROGRAM).cat

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

### find_mod_names ######

MODS    = geometry.mod globals.mod matrix.mod sparsematrix.mod 
