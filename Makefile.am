AUTOMAKE_OPTIONS = foreign
AM_INIT_AUTOMAKE = 
SUBDIRS = code_base
ACLOCAL_AMFLAGS= -I m4 
CPPFLAGS =  -I$(top_srcdir)/code_base -I$(top_srcdir)/include  -I$EIGEN_INCLUDE
CPPFLAGS +=  -O3 -funroll-loops -std=c++11
CPPFLAGS += -fopenmp -Wno-ignored-attributes -DEIGEN_NO_DEBUG -march=native

PKG_LIBS     = -lgsl -lnlopt -fopenmp 

lib_LTLIBRARIES = libDRBMD.la
libDRBMD_la_SOURCES = version.c
libDRBMD_la_LIBADD = code_base/libBMDS.la
libDRBMD_la_LIBS = -lgsl -lnlopt -fopenmp

