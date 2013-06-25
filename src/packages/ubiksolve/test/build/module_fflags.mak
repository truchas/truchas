# -*-Mode: Makefile;-*- 

ifeq ($(compiler),Absoft)

  # Absoft (whether Cygwin, Linux, or Darwin)
  MODULE_FFLAGS = -p $(libdir)
  ifneq ($(use_PGSLib),no)
    MODULE_FFLAGS += -p $(PGSLibdir)/include
  endif

else

  ifeq ($(compiler),Lahey)

    # Lahey/Fujitsu compiler (whether Cygwin or Linux)
    #
    # module files created by compilation go in the first dir in the list,
    # so put "." first
    ifneq ($(use_PGSLib),no)
      ifeq ($(os),CYGWIN)
        MODULE_FFLAGS = -mod .\;$(libdir);$(PGSLibdir)/include
      else
        MODULE_FFLAGS = --mod .:$(libdir):$(PGSLibdir)/include
      endif
    else
      ifeq ($(os),CYGWIN)
        MODULE_FFLAGS = -mod .\;$(libdir)
      else
        MODULE_FFLAGS = --mod .:$(libdir)
      endif
    endif

  else

    ifeq ($(os),SunOS)

      # Sun
      MODULE_FFLAGS = -M$(libdir)
      ifneq ($(use_PGSLib),no)
        MODULE_FFLAGS += -M$(PGSLibdir)/include
      endif

    else

      # everything else, e.g.:
      #
      # Compaq (whether Tru64Unix or Linux)
      # IBM
      # Intel (whether Cygwin or Linux)
      # NAG (whether Cygwin, Linux, SGI, DEC Alpha, Sun)
      # SGI
      MODULE_FFLAGS = -I$(libdir)
      ifneq ($(use_PGSLib),no)
        MODULE_FFLAGS += -I$(PGSLibdir)/include
      endif
    endif
  endif
endif
