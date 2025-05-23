###############################################################
#
#  First, choose a C++ compiler, and set compiler flags.
#  This is done by setting the variables CXX and CXXFLAGS.
#
###############################################################



CXX=g++
# A C++ compiler, e.g., g++, CC, xlC


CXXFLAGS=-g -O2
# Flags for the C++ compiler

CXXAUTOFLAGS= -pthread -march=native
# Flags for the C++ compiler, automatically generated by configuration script

NOCONTRACT=-ffp-contract=off -DNTL_CONTRACTION_FIXED


AR=ar
# command to make a library

ARFLAGS=ruv
# arguments for AR

RANLIB=ranlib
# set to echo if you want to disable it completely

LDFLAGS=
# libraries for linking C++ programs

LDLIBS=-lm
# libraries for linking C++ programs

CPPFLAGS=
# arguments for the C preprocessor

LIBTOOL=libtool
# libtool command -- this is now built locally

LIBTOOL_LINK_FLAGS=
# flags to add to command line when building a shared library
# mainly used to pass the argument "-no-undefined" on cygwin

LIBTOOL_LINK_LIBS=
# libraries to add to the command line when building a shared library
# mainly used to pass the argument "-lpthread", which NTL's
# configure script will do "automagically" unless user overrides
# that behavior

DEF_PREFIX=/usr/local

PREFIX=../../lib
LIBDIR=$(PREFIX)/lib
INCLUDEDIR=$(PREFIX)/include
DOCDIR=$(PREFIX)/share/doc
# where to install NTL

DESTDIR=
# added to support standard package building techniques
# that install into a "staging area"

###############################################################
#
#  Second, if you want to use GMP (the GNU Multi-Precision library),
#  define the variables GMP_OPT_INCDIR, GMP_OPT_LIBDIR, GMP_OPT_LIB below.
#  You also will have to set either NTL_GMP_LIP or NTL_GMP_HACK
#  in the config.h file.
#
#  Using GMP can lead to significant performance gains on some
#  platforms.  You can obtain GMP from http://www.swox.com/gmp.
#  Once you unpack it into a directory, just execute
#     ./configure; make
#  in that directory.
#
###############################################################


GMP_PREFIX=$(DEF_PREFIX)

GMP_INCDIR=$(GMP_PREFIX)/include
# directory containing gmp.h if using GMP 

GMP_LIBDIR=$(GMP_PREFIX)/lib
# directory containing libgmp.a if using GMP 

GMP_OPT_INCDIR=# -I$(GMP_INCDIR) # GMPI
GMP_OPT_LIBDIR=# -L$(GMP_LIBDIR) # GMPL
GMP_OPT_LIB=-lgmp # GMP
# uncomment these if using GMP


###############################################################
#
#  Third, if you want to use gf2x (a library for fast
#  multiplication over GF(2)[X]), you need to
#  define the variables GF2X_OPT_INCDIR, GF2X_OPT_LIBDIR, GF2X_OPT_LIB below.
#  You also will have to set NTL_GF2X_LIB 
#  in the config.h file.
#
###############################################################

GF2X_PREFIX=$(DEF_PREFIX)

GF2X_INCDIR=$(GF2X_PREFIX)/include
# directory containing gf2x.h if using gf2x

GF2X_LIBDIR=$(GF2X_PREFIX)/lib
# directory containing libgf2x.a

GF2X_OPT_INCDIR=# -I$(GF2X_INCDIR) # GF2X
GF2X_OPT_LIBDIR=# -L$(GF2X_LIBDIR) # GF2X
GF2X_OPT_LIB=# -lgf2x # GF2X
# uncomment these if using gf2x


###############################################################
#
# Fourth, if you do not want to run the wizard that automagically
# sets some performace related flags in config.h, set the flag below.
#
###############################################################


WIZARD=off
# Set to off if you want to bypass the wizard; otherwise, set to on.


#################################################################
#
#  That's it!  You can ignore everything else in this file!
#
#################################################################


# object files
OBJ=FFT.o FacVec.o GF2.o GF2E.o GF2EX.o GF2EXFactoring.o GF2X.o GF2X1.o  \
	GF2XFactoring.o GF2XVec.o GetTime.o GetPID.o HNF.o ctools.o LLL.o LLL_FP.o  \
	LLL_QP.o LLL_RR.o LLL_XD.o RR.o WordVector.o ZZ.o ZZVec.o ZZX.o ZZX1.o  \
	ZZXCharPoly.o ZZXFactoring.o ZZ_p.o ZZ_pE.o ZZ_pEX.o ZZ_pEXFactoring.o ZZ_pX.o  \
	ZZ_pX1.o ZZ_pXCharPoly.o ZZ_pXFactoring.o fileio.o lip.o lzz_p.o lzz_pE.o  \
	lzz_pEX.o lzz_pEXFactoring.o lzz_pX.o lzz_pX1.o lzz_pXCharPoly.o  \
	lzz_pXFactoring.o mat_GF2.o mat_GF2E.o mat_RR.o mat_ZZ.o mat_ZZ_p.o mat_ZZ_pE.o  \
	mat_lzz_p.o mat_lzz_pE.o mat_poly_ZZ.o mat_poly_ZZ_p.o mat_poly_lzz_p.o  \
	quad_float.o quad_float1.o tools.o vec_GF2.o vec_GF2E.o vec_RR.o vec_ZZ.o vec_ZZ_p.o  \
	vec_ZZ_pE.o vec_lzz_p.o vec_lzz_pE.o xdouble.o G_LLL_FP.o G_LLL_QP.o G_LLL_XD.o  \
	G_LLL_RR.o thread.o BasicThreadPool.o MatPrime.o pd_FFT.o

# library source files
SRC=FFT.cpp FacVec.cpp GF2.cpp GF2E.cpp GF2EX.cpp GF2EXFactoring.cpp GF2X.cpp  \
	GF2X1.cpp GF2XFactoring.cpp GF2XVec.cpp HNF.cpp ctools.cpp LLL.cpp LLL_FP.cpp  \
	LLL_QP.cpp LLL_RR.cpp LLL_XD.cpp RR.cpp WordVector.cpp ZZ.cpp ZZVec.cpp ZZX.cpp  \
	ZZX1.cpp ZZXCharPoly.cpp ZZXFactoring.cpp ZZ_p.cpp ZZ_pE.cpp ZZ_pEX.cpp  \
	ZZ_pEXFactoring.cpp ZZ_pX.cpp ZZ_pX1.cpp ZZ_pXCharPoly.cpp ZZ_pXFactoring.cpp  \
	fileio.cpp lip.cpp lzz_p.cpp lzz_pE.cpp lzz_pEX.cpp lzz_pEXFactoring.cpp  \
	lzz_pX.cpp lzz_pX1.cpp lzz_pXCharPoly.cpp lzz_pXFactoring.cpp mat_GF2.cpp  \
	mat_GF2E.cpp mat_RR.cpp mat_ZZ.cpp mat_ZZ_p.cpp mat_ZZ_pE.cpp mat_lzz_p.cpp  \
	mat_lzz_pE.cpp mat_poly_ZZ.cpp mat_poly_ZZ_p.cpp mat_poly_lzz_p.cpp  \
	quad_float.cpp quad_float1.cpp tools.cpp vec_GF2.cpp vec_GF2E.cpp vec_RR.cpp vec_ZZ.cpp  \
	vec_ZZ_p.cpp vec_ZZ_pE.cpp vec_lzz_p.cpp vec_lzz_pE.cpp xdouble.cpp  \
	G_LLL_FP.cpp G_LLL_QP.cpp G_LLL_XD.cpp G_LLL_RR.cpp thread.cpp  \
	BasicThreadPool.cpp MatPrime.cpp pd_FFT.cpp



# library header files
INCL=FFT.h FFT_impl.h FacVec.h GF2.h GF2E.h GF2EX.h GF2EXFactoring.h GF2X.h  \
	GF2XFactoring.h GF2XVec.h HNF.h ctools.h LLL.h RR.h WordVector.h  \
	ZZ.h ZZ_limbs.h sp_arith.h ZZVec.h ZZX.h ZZXFactoring.h ZZ_p.h ZZ_pE.h ZZ_pEX.h  \
	ZZ_pEXFactoring.h ZZ_pX.h ZZ_pXFactoring.h fileio.h lip.h lzz_p.h lzz_pE.h  \
	lzz_pEX.h lzz_pEXFactoring.h lzz_pX.h lzz_pXFactoring.h mat_GF2.h mat_GF2E.h  \
	mat_RR.h mat_ZZ.h mat_ZZ_p.h mat_ZZ_pE.h mat_lzz_p.h mat_lzz_pE.h mat_poly_ZZ.h  \
	mat_poly_ZZ_p.h mat_poly_lzz_p.h matrix.h pair.h vector.h pair_GF2EX_long.h  \
	pair_GF2X_long.h pair_ZZX_long.h pair_ZZ_pEX_long.h pair_ZZ_pX_long.h  \
	pair_lzz_pEX_long.h pair_lzz_pX_long.h quad_float.h tools.h vec_GF2.h  \
	vec_GF2E.h vec_GF2XVec.h vec_RR.h vec_ZZ.h vec_ZZVec.h vec_ZZ_p.h vec_ZZ_pE.h  \
	vec_double.h vec_long.h vec_lzz_p.h vec_lzz_pE.h vec_quad_float.h vec_vec_GF2.h  \
	vec_vec_GF2E.h vec_vec_RR.h vec_vec_ZZ.h vec_vec_ZZ_p.h vec_vec_ZZ_pE.h  \
	vec_vec_long.h vec_vec_lzz_p.h vec_vec_lzz_pE.h vec_xdouble.h xdouble.h  \
	config.h version.h new.h vec_ulong.h vec_vec_ulong.h SmartPtr.h  \
	Lazy.h LazyTable.h thread.h  BasicThreadPool.h MatPrime.h PD.h pd_FFT.h  \
	linux_s390x.h



# test data
TD=BerlekampTestIn BerlekampTestOut CanZassTestIn CanZassTestOut   \
	ZZXFacTestIn ZZXFacTestOut MoreFacTestIn LLLTestIn LLLTestOut RRTestIn RRTestOut   \
	MatrixTestIn MatrixTestOut CharPolyTestIn    \
	CharPolyTestOut QuadTestIn QuadTestOut


# test source files
TS=QuickTest.cpp ZZTest.cpp SSMulTest.cpp ZZ_pXTest.cpp lzz_pXTest.cpp BerlekampTest.cpp \
	CanZassTest.cpp ZZXFacTest.cpp  MoreFacTest.cpp \
	LLLTest.cpp subset.cpp MatrixTest.cpp mat_lzz_pTest.cpp  \
	CharPolyTest.cpp RRTest.cpp QuadTest.cpp GF2XTest.cpp GF2EXTest.cpp GF2EXGCDTest.cpp  \
	BitMatTest.cpp ZZ_pEXTest.cpp ZZ_pEXGCDTest.cpp lzz_pEXTest.cpp lzz_pEXGCDTest.cpp  \
	Timing.cpp ThreadTest.cpp ExceptionTest.cpp

# aux source to help compute crossovers
CROSS=GF2EXDivCross.cpp GF2EXGCDCross.cpp GF2EXKarCross.cpp GF2EXModCross.cpp


# scripts
SCRIPTS=MakeGetTime MakeGetPID MakeCheckFeatures ResetFeatures CopyFeatures  \
	TestScript dosify unixify RemoveProg configure DoConfig mfile cfile ppscript


# auxilliary source
MD=MakeDesc.cpp MakeDescAux.cpp newnames.cpp gen_gmp_aux.cpp gf2x_version_1_2_or_later_required.cpp
GT=GetTime0.cpp GetTime1.cpp GetTime2.cpp GetTime3.cpp GetTime4.cpp GetTime5.cpp TestGetTime.cpp
GP=GetPID1.cpp GetPID2.cpp TestGetPID.cpp
CH=CheckCompile.cpp GenConfigInfo.cpp CheckContract.cpp CheckContractAux.cpp \
	CheckThreads.cpp

AUXPROGS = TestGetTime TestGetPID CheckFeatures CheckCompile GenConfigInfo CheckContract \
	CheckThreads

FEATURES=ALIGNED_ARRAY BUILTIN_CLZL LL_TYPE SSSE3 AVX PCLMUL AVX2 FMA AVX512F \
	COPY_TRAITS1 COPY_TRAITS2 CHRONO_TIME MACOS_TIME POSIX_TIME AES_NI KMA


# documentation


DFILES=copying.txt BasicThreadPool.txt GF2.txt GF2E.txt GF2EX.txt  \
	GF2EXFactoring.txt GF2X.txt GF2XFactoring.txt GF2XVec.txt HNF.txt Lazy.txt  \
	LazyTable.txt LLL.txt RR.txt SmartPtr.txt ZZ.txt ZZ_limbs.txt ZZVec.txt ZZX.txt  \
	ZZXFactoring.txt ZZ_p.txt ZZ_pE.txt ZZ_pEX.txt ZZ_pEXFactoring.txt ZZ_pX.txt  \
	ZZ_pXFactoring.txt conversions.txt flags.txt lzz_p.txt lzz_pE.txt lzz_pEX.txt  \
	lzz_pEXFactoring.txt lzz_pX.txt lzz_pXFactoring.txt mat_GF2.txt mat_GF2E.txt  \
	mat_RR.txt mat_ZZ.txt mat_ZZ_p.txt mat_ZZ_pE.txt mat_lzz_p.txt mat_lzz_pE.txt  \
	mat_poly_ZZ.txt mat_poly_ZZ_p.txt mat_poly_lzz_p.txt matrix.txt pair.txt  \
	vector.txt quad_float.txt sedscript.txt tools.txt vec_GF2.txt vec_GF2E.txt  \
	vec_RR.txt vec_ZZ.txt vec_ZZ_p.txt vec_ZZ_pE.txt vec_lzz_p.txt vec_lzz_pE.txt  \
	xdouble.txt names.txt tour-ack.html tour-intro.html tour-time.html  \
	tour-changes.html tour-modules.html tour-unix.html tour-examples.html  \
	tour-roadmap.html tour-win.html tour-impl.html tour-struct.html tour.html  \
	tour-ex1.html tour-ex2.html tour-ex3.html tour-ex4.html tour-ex5.html  \
	tour-ex6.html tour-ex7.html arrow1.gif arrow2.gif arrow3.gif tour-gmp.html  \
	tour-gf2x.html tour-tips.html config.txt version.txt



TXFILES=GF2.txt GF2E.txt GF2EX.txt GF2EXFactoring.txt GF2X.txt  \
	GF2XFactoring.txt GF2XVec.txt HNF.txt Lazy.txt LazyTable.txt LLL.txt RR.txt  \
	SmartPtr.txt ZZ.txt ZZ_limbs.txt ZZVec.txt ZZX.txt ZZXFactoring.txt ZZ_p.txt ZZ_pE.txt  \
	ZZ_pEX.txt ZZ_pEXFactoring.txt ZZ_pX.txt ZZ_pXFactoring.txt lzz_p.txt  \
	lzz_pE.txt lzz_pEX.txt lzz_pEXFactoring.txt lzz_pX.txt lzz_pXFactoring.txt  \
	mat_GF2.txt mat_GF2E.txt mat_RR.txt mat_ZZ.txt mat_ZZ_p.txt mat_ZZ_pE.txt  \
	mat_lzz_p.txt mat_lzz_pE.txt mat_poly_ZZ.txt mat_poly_ZZ_p.txt  \
	mat_poly_lzz_p.txt matrix.txt pair.txt quad_float.txt tools.txt vec_GF2.txt  \
	vec_GF2E.txt vec_RR.txt vec_ZZ.txt vec_ZZ_p.txt vec_ZZ_pE.txt vec_lzz_p.txt  \
	vec_lzz_pE.txt vector.txt version.txt xdouble.txt BasicThreadPool.txt


HTFILES=GF2.cpp.html GF2E.cpp.html GF2EX.cpp.html GF2EXFactoring.cpp.html  \
	GF2X.cpp.html GF2XFactoring.cpp.html GF2XVec.cpp.html HNF.cpp.html  \
	Lazy.cpp.html LazyTable.cpp.html LLL.cpp.html RR.cpp.html SmartPtr.cpp.html  \
	ZZ.cpp.html ZZ_limbs.cpp.html ZZVec.cpp.html ZZX.cpp.html ZZXFactoring.cpp.html ZZ_p.cpp.html  \
	ZZ_pE.cpp.html ZZ_pEX.cpp.html ZZ_pEXFactoring.cpp.html ZZ_pX.cpp.html  \
	ZZ_pXFactoring.cpp.html lzz_p.cpp.html lzz_pE.cpp.html lzz_pEX.cpp.html  \
	lzz_pEXFactoring.cpp.html lzz_pX.cpp.html lzz_pXFactoring.cpp.html  \
	mat_GF2.cpp.html mat_GF2E.cpp.html mat_RR.cpp.html mat_ZZ.cpp.html  \
	mat_ZZ_p.cpp.html mat_ZZ_pE.cpp.html mat_lzz_p.cpp.html mat_lzz_pE.cpp.html  \
	mat_poly_ZZ.cpp.html mat_poly_ZZ_p.cpp.html mat_poly_lzz_p.cpp.html  \
	matrix.cpp.html pair.cpp.html quad_float.cpp.html tools.cpp.html  \
	vec_GF2.cpp.html vec_GF2E.cpp.html vec_RR.cpp.html vec_ZZ.cpp.html  \
	vec_ZZ_p.cpp.html vec_ZZ_pE.cpp.html vec_lzz_p.cpp.html vec_lzz_pE.cpp.html  \
	vector.cpp.html version.cpp.html xdouble.cpp.html BasicThreadPool.cpp.html



DOC = $(DFILES) $(HTFILES) TFT-time.jpg zmulrat.jpg flintrat.jpg


# test program executables
PROGS=QuickTest ZZTest SSMulTest ZZ_pXTest lzz_pXTest BerlekampTest CanZassTest \
	ZZXFacTest MoreFacTest LLLTest  \
	BitMatTest MatrixTest mat_lzz_pTest CharPolyTest RRTest QuadTest GF2XTest  \
	GF2EXTest GF2EXGCDTest subset ZZ_pEXTest ZZ_pEXGCDTest lzz_pEXTest lzz_pEXGCDTest \
	Timing ThreadTest

# things to save to a tar file
SFILES=$(SRC) $(SCRIPTS) $(MD) $(GT) $(GP) $(CH) $(TS) $(TD) $(CROSS) mach_desc.win   \
	Poly1TimeTest.cpp Poly2TimeTest.cpp Poly3TimeTest.cpp GF2XTimeTest.cpp    \
	InitSettings.cpp DispSettings.cpp  WizardAux Wizard 


#################################################################
#
#  Rules for compiling the library
#
#################################################################


NTL_INCLUDE = -I../include -I.
# NTL needs this to find its include files

COMPILE = $(CXX) $(NTL_INCLUDE) $(CPPFLAGS) $(CXXFLAGS) $(CXXAUTOFLAGS)  -c

LINK = $(CXX) $(NTL_INCLUDE) $(CPPFLAGS) $(CXXFLAGS) $(CXXAUTOFLAGS) $(LDFLAGS)



# 'make' or 'make all' does a complete make, including additional
# setup not done in configure.

# The file setup-phase is removed by the configure script
# when it terminates successfully.

# The file need-to-run-configure is removed by the configure script
# before making any changes to makefile/config.h and is recreated
# when it terminates successfully.

all: setup-phase
	$(MAKE) ntl.a

setup-phase: need-to-run-configure
	$(MAKE) clobber
	$(MAKE) setup1
	$(MAKE) setup2
	$(MAKE) setup3
	$(MAKE) setup4
	touch setup-phase


# setup1 generates the file ../incluse/NTL/mach_desc.h

setup1: 
	$(COMPILE) MakeDescAux.cpp
	$(LINK) -o MakeDesc MakeDesc.cpp MakeDescAux.o $(LDLIBS)
	./MakeDesc
	mv mach_desc.h ../include/NTL/mach_desc.h


# setup2 does some dynamic checks for GetTime, GetPID, and other features

setup2:
	echo "*** CheckFeatures log ***" > CheckFeatures.log
	sh MakeGetTime "$(LINK)" "$(LDLIBS)"
	sh MakeGetPID  "$(LINK)" "$(LDLIBS)"
	sh MakeCheckFeatures "$(FEATURES)" "$(LINK)" "$(LDLIBS)"

# NOTE: to add a feature XXX:
#   * add a program CheckXXX.cpp which returns 0 if XXX works, -1 otherwise
#   * add XXX to the FEATURES variable

# setup3 generates the file ../include/NTL/gmp_aux.h
# The file ../include/NTL/gmp_aux.h is included in ../include/NTL/lip.h
# when NTL_GMP_LIP is set.
# When this flag is not set, an empty files produced.
# This also checks that the right version of gf2x library.

setup3:
	$(LINK) $(GMP_OPT_INCDIR) -o gen_gmp_aux gen_gmp_aux.cpp $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(LDLIBS)
	./gen_gmp_aux > ../include/NTL/gmp_aux.h 
	$(LINK) $(GF2X_OPT_INCDIR) -o gf2x_version_1_2_or_later_required gf2x_version_1_2_or_later_required.cpp $(GF2X_OPT_LIBDIR) $(GF2X_OPT_LIB) $(LDLIBS) 

# setup4 runs the wizard

setup4:
	sh Wizard $(WIZARD) "$(MAKE)" "$(FEATURES)"


ntl.a:	$(OBJ) 
	$(AR) $(ARFLAGS) ntl.a $(OBJ) #LSTAT
	- $(RANLIB) ntl.a #LSTAT
# 	$(LIBTOOL) --tag=CXX --mode=link $(LINK) $(LIBTOOL_LINK_FLAGS) -o libntl.la $(OBJ:.o=.lo) $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(GF2X_OPT_LIBDIR) $(GF2X_OPT_LIB) $(LDLIBS) $(LIBTOOL_LINK_LIBS) -rpath $(LIBDIR) -version-info `cat VERSION_INFO` #LSHAR

LCOMP= #LSTAT
# LCOMP=$(LIBTOOL) --tag=CXX --mode=compile #LSHAR

lip.o:	lip.cpp 
	$(LCOMP) $(COMPILE) $(GMP_OPT_INCDIR) lip.cpp

GF2X.o:	GF2X.cpp
	$(LCOMP) $(COMPILE) $(GF2X_OPT_INCDIR) GF2X.cpp

quad_float.o:	quad_float.cpp
	$(LCOMP) $(COMPILE) $(NOCONTRACT) quad_float.cpp

pd_FFT.o: pd_FFT.cpp
	$(LCOMP) $(COMPILE) $(NOCONTRACT) pd_FFT.cpp

CheckCompile: CheckCompile.cpp
	$(LINK) -o  CheckCompile CheckCompile.cpp $(LDLIBS)

GenConfigInfo: GenConfigInfo.cpp
	$(LINK) -o  GenConfigInfo GenConfigInfo.cpp $(LDLIBS)

CheckContract: CheckContract.cpp CheckContractAux.cpp
	$(LINK) $(NOCONTRACT) -o  CheckContract CheckContract.cpp CheckContractAux.cpp $(LDLIBS)

CheckThreads: CheckThreads.cpp
	$(LINK) -o  CheckThreads CheckThreads.cpp $(LDLIBS)


.cpp.o: 
	$(LCOMP) $(COMPILE) $<

.cpp: 
	$(LINK) -o $@ $< ntl.a $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(GF2X_OPT_LIBDIR) $(GF2X_OPT_LIB) $(LDLIBS) #LSTAT
# 	$(LIBTOOL) --tag=CXX --mode=link $(LINK) -o $@ $< libntl.la #LSHAR

#################################################################
#
# Rule for running tests
# make check runs a series of tests
#
#################################################################

check:
	sh RemoveProg $(PROGS)
	$(MAKE) QuickTest
	./QuickTest
	sh RemoveProg QuickTest
	sh TestScript "$(MAKE)"

#################################################################
#
# Rule for installing
# make install just does a simple copy of the include file
#   and library.  The -p option is used to preserve file attributes.
#   This helps avoid some problems (especially when copying ntl.a).
#   Also, an attempt is made to make everything that is
#   installed readable by everyone.
#   
# make uninstall removes these files
#   
#################################################################




install:
	mkdir -p -m 755 $(DESTDIR)$(INCLUDEDIR)
	rm -rf $(DESTDIR)$(INCLUDEDIR)/NTL
	mkdir -m 755 $(DESTDIR)$(INCLUDEDIR)/NTL
	cp -p ../include/NTL/*.h $(DESTDIR)$(INCLUDEDIR)/NTL
	- chmod -R a+r $(DESTDIR)$(INCLUDEDIR)/NTL
	mkdir -p -m 755 $(DESTDIR)$(DOCDIR)
	rm -rf $(DESTDIR)$(DOCDIR)/NTL
	mkdir -m 755 $(DESTDIR)$(DOCDIR)/NTL
	cp -p ../doc/*.txt $(DESTDIR)$(DOCDIR)/NTL
	cp -p ../doc/*.html $(DESTDIR)$(DOCDIR)/NTL
	cp -p ../doc/*.gif $(DESTDIR)$(DOCDIR)/NTL
	- chmod -R a+r $(DESTDIR)$(DOCDIR)/NTL
	mkdir -p -m 755 $(DESTDIR)$(LIBDIR)
	cp -p ntl.a $(DESTDIR)$(LIBDIR)/libntl.a #LSTAT
	- chmod a+r $(DESTDIR)$(LIBDIR)/libntl.a #LSTAT
# 	$(LIBTOOL) --mode=install cp -p libntl.la $(DESTDIR)$(LIBDIR) #LSHAR


uninstall:
	rm -f $(LIBDIR)/libntl.a #LSTAT
# 	$(LIBTOOL) --mode=uninstall rm -f $(LIBDIR)/libntl.la #LSHAR
	rm -rf $(INCLUDEDIR)/NTL
	rm -rf $(DOCDIR)/NTL

#################################################################
#
#  Rules for cleaning up
#
#  make clobber removes *everything* created by make,
#    but it does not restore config.h to its default.
#
#  make clean tidies up a bit
#
#################################################################

clobber:	
	rm -f ntl.a mach_desc.h ../include/NTL/mach_desc.h  GetTime.cpp GetPID.cpp
	sh ResetFeatures '..' "$(FEATURES)" 
	rm -f ../include/NTL/gmp_aux.h
	sh RemoveProg $(PROGS) MakeDesc $(AUXPROGS) gen_gmp_aux gf2x_version_1_2_or_later_required
	rm -f *.o
	rm -rf small
	rm -f cfileout mfileout
	rm -rf .libs *.lo libntl.la
	rm -f setup-phase

clean:
	sh RemoveProg $(PROGS) MakeDesc $(AUXPROGS) gen_gmp_aux gf2x_version_1_2_or_later_required
	rm -f *.o 
	rm -rf small
# 	- $(LIBTOOL) --mode=clean rm -f libntl.la *.lo #LSHAR

wclean:
	rm -f *.o 

#################################################################
#
#  Rules for making tar and zip files
# 
#  make libtool-gen-origin generates the directory 
#     libtool-origin used to include in the distribution
#     - this only needs to be run very occasionally, to keep
#       libtool relatively up-to-date
#     - it must be run on a machine with autotools
#
#  make ppdoc creates pretty-printed versions of some documentation
#     - run before make package or make winpack
#
#  make package creates a tar.gz file suitable for Unix
#
#  make winpack creates a zip file suitable for Windows
#
#################################################################

libtool-gen-origin:
	rm -rf libtool-origin && \
	   cp -R libtool-seed libtool-origin && \
	   cd libtool-origin && autoreconf -fiv && rm -rf autom4te.cache

ppdoc:
	sh ppscript "$(TXFILES)"

ppclean:
	rm -f ../doc/*.cpp


package:
	sh unixify "$(SFILES) DIRNAME WINDIR VERSION_INFO NOTES" "$(INCL)" "$(DOC)" "$(FEATURES)"
	rm -rf `cat DIRNAME`
	rm -f `cat DIRNAME`.tar
	rm -f `cat DIRNAME`.tar.gz
	mv unix `cat DIRNAME`
	chmod -R a+rX `cat DIRNAME`
	tar -cvf `cat DIRNAME`.tar `cat DIRNAME`
	gzip `cat DIRNAME`.tar
	rm -rf `cat DIRNAME`

winpack:
	./configure --nowrite NTL_GMP_LIP=off NTL_TLS_HACK=off
	sh dosify "$(SRC)" "$(INCL)" "$(DOC)" "$(TS)" "$(TD)" "$(FEATURES)"
	rm -rf `cat WINDIR`
	rm -f `cat WINDIR`.zip 
	mv dos `cat WINDIR`
	chmod -R a+rX `cat WINDIR`
	find ./`cat WINDIR` '!' '(' -name '*.gif' -o -name '*.jpg' ')'  -print | zip -l `cat WINDIR` -@
	find ./`cat WINDIR` -name '*.gif' -o -name '*.jpg' -print | zip -u `cat WINDIR` -@
	rm -rf `cat WINDIR`


######################################################################
#
#  config wizard related stuff
#
######################################################################

WOBJ=FFT.o GetTime.o GetPID.o ctools.o ZZ.o ZZVec.o ZZ_p.o ZZ_pX.o ZZ_pX1.o  \
	lip.o tools.o vec_ZZ.o vec_ZZ_p.o GF2.o WordVector.o vec_GF2.o GF2X.o GF2X1.o  \
	thread.o BasicThreadPool.o fileio.o


# wntl.a: LCOMP= #LSHAR
wntl.a:	$(WOBJ) 
	$(AR) $(ARFLAGS) wntl.a $(WOBJ)
	- $(RANLIB) wntl.a


Poly1TimeTest:
	$(LINK) -o Poly1TimeTest Poly1TimeTest.cpp wntl.a  $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(LDLIBS)
Poly2TimeTest:
	$(LINK) -o Poly2TimeTest Poly2TimeTest.cpp wntl.a  $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(LDLIBS)
Poly3TimeTest:
	$(LINK) -o Poly3TimeTest Poly3TimeTest.cpp wntl.a  $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(LDLIBS)


GF2XTimeTest:
	$(LINK) -o GF2XTimeTest GF2XTimeTest.cpp wntl.a  $(GMP_OPT_LIBDIR) $(GMP_OPT_LIB) $(LDLIBS)

InitSettings:
	$(LINK) -o InitSettings InitSettings.cpp  $(LDLIBS)


DispSettings:
	$(LINK) -o DispSettings DispSettings.cpp $(LDLIBS)

