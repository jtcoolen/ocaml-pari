dir=examples
file=$objdir/Makefile.examples
lnfile=Makefile

echo Extracting $file
rm -f $file $dir/$lnfile

RUNPTH=
if test -n "$runpathprefix"; then
  RUNPTH=`config/ldflags "$LDneedsWl" "$runpathprefix $runpath"`
fi

cat > $file << EOT
# Generic Makefile for PARI programs -- $pretty
#
#  This file was created by Configure. Any change made to it will be
#  lost when Configure is run.
#
# make all will create
#  extgcd (linked with libpari)
#  libextgcd.so (to be used by "install" under GP)
#
# Under GP: install("extgcd", "GG&&", "f", "./libextgcd.so") enables
# you to subsequently use f(a,b) to call extgcd(a,b)

# change this TARGET to compile your own programs
TARGET = extgcd
SHELL  = $make_sh

DBGFLAGS   = $DBGFLAGS
CFLAGS     = $OPTFLAGS
EXTRACFLAGS=
#CFLAGS    = \$(DBGFLAGS)

CC         = $CC
CPPFLAGS   = -I. -I$includedir
LD         = $LD
LDFLAGS    = $LDFLAGS
MODLD      = $MODLD
MODLDFLAGS = $MODLDFLAGS
EXTRAMODLDFLAGS = $EXTRAMODLDFLAGS
EXTRALIBS  =

RUNPTH     = $RUNPTH
DLCFLAGS   = $DLCFLAGS
LIBS       = $LIBS -L$libdir -lpari

RM = rm -f

OBJS = \$(TARGET).o
DYN = lib\$(TARGET).$DLSUFFIX
ALL = \$(TARGET) \$(DYN)

all: \$(ALL)

dynlib: \$(DYN)

\$(DYN): \$(OBJS)
	\$(MODLD) -o \$@ \$(MODLDFLAGS) \$(EXTRACFLAGS) \$(OBJS) \$(EXTRAMODLDFLAGS)

\$(TARGET): \$(OBJS)
	\$(LD) -o \$@ \$(LDFLAGS) \$(EXTRACFLAGS) \$< \$(EXTRALIBS) \$(RUNPTH) \$(LIBS)

%.o: %.c
	\$(CC) -c \$(CFLAGS) \$(EXTRACFLAGS) \$(CPPFLAGS) \$(DLCFLAGS) \$<
clean:
	-\$(RM) *.o \$(ALL)
EOT

( cd $dir ; $ln_s ../$file $lnfile )
