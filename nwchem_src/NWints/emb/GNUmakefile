# $Id: GNUmakefile 22012 2012-02-12 16:52:25Z jhammond $

    OBJ          = 
    OBJ_OPTIMIZE = embini.o embint.o

  USES_BLAS = embini.o embint.o

  LIB_TARGETS = embini.o embint.o
  LIBRARY = libnwints.a

  HEADERS= embcom.fh

include ../../config/makefile.h
include ../../config/makelib.h

embini.o:	embini.F	embcom.fh
embint.o:   embint.F	embcom.fh
