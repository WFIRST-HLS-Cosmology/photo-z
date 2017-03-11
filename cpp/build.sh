#! /bin/bash

# temporary build script for use until a Makefile is created


# for looking at assembly: -S -fverbose-asm -masm=intel 
# -Wa,-alh will write assembly and source to std out.  

OPTS="-std=c++11 -O3 -march=native -ftree-vectorize -fdata-sections -fno-rtti -Wno-unused-result"
INCLUDES="-I ../3rd-party/include/"
FMTPATH="../3rd-party/include/fmt"
LIBS="-L ./"

utils="g++ $OPTS -c basic_utils.cpp $FMTPATH/format.cc $INCLUDES" 
archive="ar rvs libutils.a format.o basic_utils.o"
makesim="g++ $OPTS -o make-grid make_grid.cpp sed.cpp grid_tools.cpp libutils.a $INCLUDES $LIBS -lgsl -lgslcblas -lm -fopenmp"
photz="g++ $OPTS -o photo-z photo_z.cpp fitting_tools.cpp libutils.a $INCLUDES $LIBS -lm -fopenmp"
fitcat="g++ $OPTS -o fitcat fit_catalog.cpp grid_tools.cpp sed.cpp libutils.a $INCLUDES $LIBS -lgsl -lgslcblas -lm -fopenmp"

#g++ $OPTS -fexceptions -S -fverbose-asm -masm=intel fitting_tools.cpp $INCLUDES -lm -fopenmp

if [ -e libutils.a ]; then 
    if [ libutils.a -ot basic_utils.cpp ] || [ libutils.a -ot basic_utils.hpp ] ; then
        $($utils) && $archive
    fi
else
    $($utils) && $archive
fi

$($makesim) & 
$($photz) & 
$($fitcat) 

wait

strip photo-z make-grid fitcat

if [ "$STRATOS_HOME" != "" ]; then

    echo "Copying executables to $STRATOS_HOME"
    sudo cp photo-z $STRATOS_HOME/bin/
    sudo cp make-grid $STRATOS_HOME/bin/
    sudo cp fitcat $STRATOS_HOME/bin/
fi
