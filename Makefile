#----------------------------------------------------------------------------------------------------
# simple Makefile for T-Digest, relies on cmake to do the actual build.  
# Use CMAKE_LIBRARY_OPTIONS or CMAKE_FULL_OPTIONS argument to this Makefile to pass options to cmake.
#----------------------------------------------------------------------------------------------------

ifndef CMAKE_LIBRARY_OPTIONS
	CMAKE_LIBRARY_OPTIONS=\
		-DBUILD_SHARED=ON \
		-DENABLE_FRAME_POINTER=ON \
		-DENABLE_FRAME_POINTER=ON \
		-DBUILD_TESTS=OFF \
		-DBUILD_EXAMPLES=OFF
endif

ifndef CMAKE_FULL_OPTIONS
	CMAKE_FULL_OPTIONS=\
		-DBUILD_SHARED=ON \
		-DENABLE_FRAME_POINTER=ON \
		-DENABLE_FRAME_POINTER=ON \
		-DBUILD_TESTS=ON \
		-DBUILD_EXAMPLES=ON
endif

default: full

# just build the static and dynamic libraries. Do not build tests or benchmarks
library:
	( cd build ; cmake $(CMAKE_LIBRARY_OPTIONS) .. ; $(MAKE) )

# build all
full:
	( cd build ; cmake $(CMAKE_FULL_OPTIONS) .. ; $(MAKE) )

clean: distclean

distclean:
	rm -rf build/* 
