CC            = pgcc

CFLAGS        = -fast

FC            = pgf90

FFLAGS        = -fast -Mscalarsse -Mvect=sse,cachesize:262144 -Mflushz -pc 64 -Mcache_align -O4
# Per compilare codici > 2GB
#FFLAGS        = -O2 -mcmodel=medium
FFLAGS        = -fast

DEST	      = .

EXTHDRS	      =

HDRS	      =

# Per compilare codici > 2GB (e usare lapack customizzate)
#LDFLAGS	      = -L. -L/usr/local/lib -mcmodel=medium

LDFLAGS	      = -L. -L/usr/local/lib \
                -L/usr/local/pgi/linux86-64/7.1/lib 

LIBS	      = /usr/lib/liblapack.a /usr/lib/libblas.a

LINKER	      = $(FC)

MAKEFILE      = Makefile

OBJS	      = nodipiante.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = nodipiante.f

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

debug:;         @make "FFLAGS=-g" $(PROGRAM)

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

