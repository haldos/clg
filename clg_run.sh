#!/bin/bash

EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {arg}"
  exit 1
fi

pathfilename1=${1%.*}
fileext1=${1##*.}
pathfilename2=${2%.*}
fileext2=${2##*.}
pathname=${pathfilename1%/*}

# PARAMETERS
ALPHA=100
NUMIT=50
SIGMA=0.5
RHO=30

# GRAYSCALE CONVERSION
BW1="$pathfilename1.bw.$fileext1"
bin/test_rgb2gray $1 $BW1
echo RGB image $1 converted to grayscale image "$pathfilename1.bw.$fileext1".
BW2="$pathfilename2.bw.$fileext2"
bin/test_rgb2gray $2 $BW2
echo RGB image $2 converted to grayscale image "$pathfilename2.bw.$fileext2".
echo Grayscale conversion done.

# CLG RUN
TIFF="$pathname/flow.tiff"
bin/clg_of $BW1 $BW2 $ALPHA $RHO $SIGMA $NUMIT $TIFF
echo CLG calculation done, output saved in $TIFF.

# TIFF SWAPPING
bin/test_swaptiff $TIFF
echo TIFF file $TIFF swapped.

# VIEWFLOW
PNG="$pathname/flow.png"
../../imscript64/bin/viewflow -1 $TIFF $PNG
echo PNG visualization generated, saved in $PNG.

# REMOVE BW & TIFF IMAGES
rm $BW1 $BW2 $TIFF
echo Removing $BW1 $BW2 $TIFF... Removed.
