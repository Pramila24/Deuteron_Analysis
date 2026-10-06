# make file for F2PY and installation into the utilities 
#
# default procedure for compiling a .f file
#
#suffix rule necessary for compiling with absoft
.f.so :
	$(F2PY) $(F2PYFLAGS) $<

HERE          = ./
PYTHON        = python3

CROSEC        = ./

CERNLIB       =

DEST          = ./build/


PRINT	      = 

F2PY          = $(PYTHON) -m numpy.f2py
F2PYFLAGS     = -c --build-dir ./build/ -m

PROGRAM	      = Laget_Xsec_fp_1
EXT_SUFFIX    = $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
PROGRAM_EXT   = $(PROGRAM)$(EXT_SUFFIX)

SRCS          = Laget_Xsec_fp_1.f

all:		$(PROGRAM_EXT)

$(PROGRAM_EXT): $(SRCS)
		$(F2PY) $(F2PYFLAGS) $(PROGRAM)  $(SRCS)

clean:	
		rm -f $(PROGRAM)*.so
		rm -rf ./build
		rm -f *.o








