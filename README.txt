IPOL 2D CLG Optical Flow

This program implements the combined local-global (CLG) algorithm for optical 
flow computation in 2D greyscale images, introduced by Bruhn et al. (2002). 
The code is supplied in addition to the IPOL documentation material, available 
online at http://www.ipol.im.

Program written by Jorge Jara <jjara@dcc.uchile.cl> and Jose Delpiano 
<jdelpian@uandes.cl>. Distributed under the terms of the LGPLv3 license.


== Usage ==


From command line:

clg_of <image1> <image2> <alpha> <rho> <sigma> <iterations> <flow_output_image>

<image1>, <image2> specify the input time frames (images) for the optical flow 
                   computation

<alpha>            specifies the value of the global regularization coefficient.
                   Use a value greater than 0.

<rho>              specifies the value of the local spatio-temporal derivatives 
                   regularization coefficient in the form of a gaussian 
                   smoothing with standard deviation rho (greater than 0). Use 
                   0 to avoid this smoothing.

<sigma>            specifies standard deviation value (greater than 0) for a 
                   gaussian smoothing applied to the input images. Use 0 to 
                   avoid smoothing.


== Compiling ==

Instructions are included below for compiling on Linux sytems with GCC, on 
Windows with MSVC.


== Compiling (Linux) ==

To compile this software under Linux, first install the development files for
libjpeg, libpng, and libtiff.  On Ubuntu and other Debian-based systems, enter
the following into a terminal:
    sudo apt-get install build-essential libjpeg8-dev libpng-dev libtiff-dev
On Redhat, Fedora, and CentOS, use
    sudo yum install make gcc libjpeg-turbo-devel libpng-devel libtiff-devel

--- Building with PNG, JPEG, and/or TIFF support ---

To use the program with PNG, JPEG, and/or TIFF images, the 
following libraries are needed.

    For PNG:    libpng and zlib
    For JPEG:   libjpeg 
    For TIFF:   libtiff

These libraries can be obtained at 
    
    http://www.libpng.org/pub/png/libpng.html
    http://www.zlib.net/
    http://www.ijg.org/
    http://www.remotesensing.org/libtiff/


== Compiling (Windows with MSVC) ==

The express version of the Microsoft Visual C++ (MSVC) compiler can be 
obtained for free at

    http://www.microsoft.com/visualstudio/en-us/products/2010-editions/express

To use the program with PNG, JPEG, and/or TIFF images, the 
following libraries are needed.

    For PNG:    libpng and zlib
    For JPEG:   libjpeg 
    For TIFF:   libtiff

You will have to specify the location of the corresponding libraries when 
compiling.

== Acknowledgements ==

This material is based upon work supported by Comisión Nacional de 
Investigación Científica y Tecnológica (CONICYT), Chile, project No. 1090246. 
Any opinions, findings, and conclusions or recommendations expressed in this 
material are those of the author(s) and do not necessarily reflect the views 
of the CONICYT. Jorge Jara and Jose Delpiano hold PhD scholarships granted by 
CONICYT.
