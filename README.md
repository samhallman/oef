# Oriented Edge Forests

This repository contains the source code for the method described in

"Oriented Edge Forests for Boundary Detection".
Sam Hallman, Charless Fowlkes. CVPR, June 2015

The system is implemented in MATLAB. On a 480-by-320 image, the detector should
run in ~2 seconds on an 8-core machine. Development was done on Linux and
pre-compiled MEX binaries for Linux are included.


## Installation

To use this software, you need to have Piotr Dollar's very useful
[Image & Video Matlab Toolbox] installed.

You can download a pre-trained model at
http://www.ics.uci.edu/~shallman/oef/modelCvpr.mat.
The file is 98 MB, but swells to 1.1 GB when loaded into memory.
To train a model yourself, you'll need to download the [BSDS500] dataset.


## Usage

See `demo.m` for usage examples.

To train a reasonably good detector quickly,

    % requires ~5GB of RAM and <4 min/tree
    model = train('nPos',5e5, 'nNeg',5e5, 'nTrees',8, ...
      'useParfor',1, 'calibrate',0, 'bsdsDir','/path/to/bsds/');

To train the model from the CVPR paper, just use the default settings:

    % requires ~19GB of RAM and ~15 min/tree
    model = train('bsdsDir','/path/to/bsds/');

This trains 24 trees by default, because that is originally how I derived the
numbers shown in the paper. But 24 trees is probably overkill, and I would bet
that you'd get the same results with 12 trees.

### Acknowledgements

Many files were built on top of files from the [Sketch Tokens] and [Structured
Forest] packages. I also make use of the edge linking files from Peter Kovesi's
[MATLAB and Octave Functions for Computer Vision and Image Processing] page.


[Sketch Tokens]:https://github.com/joelimlimit/SketchTokens
[Structured Forest]:https://github.com/pdollar/edges
[BSDS500]:http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html#bsds500
[Image & Video Matlab Toolbox]:http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html
[MATLAB and Octave Functions for Computer Vision and Image Processing]:http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/index.html
