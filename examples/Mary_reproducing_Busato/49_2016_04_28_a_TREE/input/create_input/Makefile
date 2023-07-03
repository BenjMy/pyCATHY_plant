CC            = gcc

CFLAGS        = -fast

FC            = gfortran

# Per compilare codici > 2GB
FFLAGS        = -O5 -mcmodel=medium

DEST	      = .

EXTHDRS	      =

HDRS	      =

# Per compilare codici > 2GB (e usare lapack customizzate)
LDFLAGS	      = -L. -L/usr/local/lib -mcmodel=medium

LIBS	      = -llapack -lblas

LINKER	      = $(FC)

MAKEFILE      = Makefile

OBJS	      = meteoINPUT.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = out_meteo.exe

SRCS	      = meteoINPUT.f

.SUFFIXES: .prj .f90

.f.prj:
	ftnchek -project -declare -noverbose $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

PRJS= $(SRCS:jo=.prj)

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)

check: $(PRJS)
	ftnchek -noverbose -declare $(PRJS) -project -noextern -library > $(PROGRAM).ftn

profile:;       @make "FFLAGS=-pg" "CFLAGS=-pg" "LDFLAGS=-pg" $(PROGRAM)

debug:;         @make "FFLAGS=-g -C" $(PROGRAM)

clean:;		rm -f $(OBJS) $(PROGRAM).cat $(PROGRAM).ftn

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

print:;		$(PRINT) $(HDRS) $(SRCS)

cat:;		$(CAT) $(HDRS) $(SRCS) > $(PROGRAM).cat

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

