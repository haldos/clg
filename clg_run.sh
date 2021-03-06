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
rm flow_u.png flow_v.png

# GENERATE TILED IMAGE
TRUEFLOW="$pathname/t.png"
ORIGINAL1="$pathname/1.png"
ORIGINAL2="$pathname/2.png"
TRUEFLOW_LABELED="$pathname/3.png"
FLOW_LABELED="$pathname/4.png"
TILED="$pathname/output.png"
convert $1 -gravity North -background White -splice 0x18 -annotate +0+2 'Input image 1' $ORIGINAL1
convert $2 -gravity North -background White -splice 0x18 -annotate +0+2 'Input image 2' $ORIGINAL2
convert $TRUEFLOW -gravity North -background White -splice 0x18 -annotate +0+2 'True flow' $TRUEFLOW_LABELED
convert $PNG -gravity North -background White -splice 0x18 -annotate +0+2 'CLG computed flow' $FLOW_LABELED

montage $ORIGINAL1 $ORIGINAL2 $TRUEFLOW_LABELED $FLOW_LABELED -tile 2x2 -geometry 256x274+0+0 $TILED
rm $ORIGINAL1 $ORIGINAL2 $TRUEFLOW_LABELED $FLOW_LABELED
display $TILED &

# Haldo Spontón, 2012.
# IIE - FING - UDELAR.
# haldos@fing.edu.uy

