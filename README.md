# deformer
A tool to deform images

## Sample videos
[A mosaic heated up (pixels offset randomly) and let to cool back down to straighten itself.](https://youtu.be/GvXnWSNPdg8)

[A checkerboard spun by 90-degrees in response to a subset of pixels being rotated and fixed in place.](https://youtu.be/8SkTvcZNEnU)

## Basics
This is a tool to deform images. Parts of an image can be offset, rotated and scaled, they can be fixed in place or left free to relax, and their elastic stiffness can be controlled. Mask images can be used to select the pixels to be manipulated. Individual pixels can also be manipulated.

Unfortunately there is no interactive graphical interface as I have not yet familiarized myself how to do such things using C++. I am more interested in the underlying algorithms anyway. However, there are plenty of ready-made demonstrations the user can modify to their liking; see the 'Demos' section below. Use the provided pngreader and pngwriter tools to convert .png images to .rgb input files and .rgb output files back to .png images.

Deformations are based on an underlying square grid of nodes connected by spring forces. More precisely the nodes are in a multiscale grid to provide more resolution where it is needed and to minimize the computing effort elsewhere. The elastic system is relaxed using simple gradient descent. The code is written in C++ and exploits OpenMP parallelization where reasonable, but it is not particularly high performance.

## Demos
...

## Relaxation options
...

## Multigrid
...

## Visualization
...
