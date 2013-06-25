# -*-Mode: Makefile;-*- 

# required to be set:
#   LIB
#   LIBNAME
#   libdir

LDFLAGS = $(empty)

ifeq ($(os),CYGWIN)
  ifeq ($(compiler),Intel)
    LDFLAGS += /link
  endif
  ifeq ($(compiler),Lahey)
    LDFLAGS += -Lib 
  endif
  LDFLAGS += $(LIB).lib
else
  LDFLAGS += $(STD_FFLAGS)
  LDFLAGS += -L$(libdir) -l$(LIBNAME)
  ifeq ($(DEBUG),yes)
    LDFLAGS += -g
  endif
endif

# Intel compiler
ifeq ($(compiler),Intel)
  ifeq ($(os),CYGWIN)
    ifeq ($(compiler_ver),-8.0)
      LDFLAGS += /LIBPATH:C:/Program\ Files/Intel/compiler80/IA32/lib
      LDFLAGS += /LIBPATH:C:/Program\ Files/Microsoft\ Visual\ Studio\ .NET/Vc7/lib
      LDFLAGS += /LIBPATH:C:/Program\ Files/Microsoft\ Visual\ Studio\ .NET/Vc7/PlatformSDK/lib
    else
      LDFLAGS += /LIBPATH:C:/Program\ Files/Intel/Compiler70/IA32/lib
      LDFLAGS += /LIBPATH:C:/Program\ Files/Microsoft\ Visual\ Studio/VC98/Lib
    endif
    # desired
    #  LDFLAGS = /LIBPATH:$(INTELDIR)/lib /LIBPATH:$(MSDIR)/vc98/Lib
    ifeq ($(DEBUG),yes)
      LDFLAGS += /DEBUG
    endif
  endif
endif

# Absoft compiler
ifeq ($(compiler),Absoft)
  ifeq ($(os),CYGWIN)
    LDFLAGS += $(ABSOFT)/lib/unix.lib
    LDFLAGS += -aliases:$(ABSOFT)/lib/unicode.als
    ifeq ($(DEBUG),yes)
      LDFLAGS += -g
    endif
  else
    LDFLAGS += -lU77
  endif
endif

ifeq ($(os),IRIX)
  ifneq ($(DEBUG),yes)
    LDFLAGS += -IPA
  endif
endif

ifeq ($(os),UNICOS)
  LDFLAGS += -l perf -l sci
endif
