# deformer
A tool to deform images

## Sample videos
[A mosaic heated up (pixels offset randomly) and let to cool back down to straighten itself.](https://youtu.be/GvXnWSNPdg8)

[A checkerboard spun by 90-degrees in response to a subset of pixels being rotated and fixed in place.](https://youtu.be/8SkTvcZNEnU)

## Basics
This is a tool to deform images. Parts of an image can be offset, rotated and scaled, they can be fixed in place or left free to relax, and their elastic stiffness can be controlled. Mask images can be used to select the pixels to be manipulated. Individual pixels can also be manipulated.

Unfortunately there is no interactive graphical interface as I have not yet familiarized myself how to do such things using C++. I am more interested in the underlying algorithms anyway. However, there are plenty of ready-made demonstrations the user can modify to their liking; see the 'Demos' section below.

Deformations are based on an underlying square grid of nodes connected by spring forces. More precisely, the nodes are in a multiscale grid to provide more resolution where it is needed and to minimize the computing effort elsewhere. The elastic system is relaxed using simple gradient descent. The code is written in C++ and exploits OpenMP parallelization where reasonable, but it is not particularly high performance.

## Demos
deformer/source/main.cpp contains plenty of ready-made demonstrations of how to use the code and that the user can modify to their liking. Use the provided pngreader and pngwriter tools (under deformer/source) to convert .png images to .rgb input files and .rgb output files back to .png images; instructions are given at the beginning of the code files. Some input files are provided under deformer/demos; otherwise links to recommended input images are given in deformer/source/main.cpp.

## Relaxation options
The basic gradient descent relaxation for the elastic systems comes with two modifiers: to eliminate drift and to limit step size. Elimination of drift means that any rigid body motions and rotations resulting from a relaxation step are eliminated. Note that Newton's third law is not exactly satisfied by the implementation due to simpler code and higher performance. Limiting the step size means that the updates to the nodes' positions are capped. This can help with stability issues in some instances. More precisely, the magnitude of the position updates is modulated monotonically by a hyperbolic tangent function. The coefficient allows to scale the maximum update magnitude with respect to the mean update magnitude.

## Forces on nodes
The following diagram depicting a 3-by-3 grid of nodes helps demonstrate the forces on the nodes. Note that each node has four neighbors (if not at the edge).

<img src="https://github.com/petenez/deformer/blob/main/springs.png" width="400">

For a node <img src="https://latex.codecogs.com/svg.image?n">, the force on it is given by

<img src="https://latex.codecogs.com/svg.image?F_n=\sum_{m}\varsigma_{nm}\left(\left\lbrace\begin{matrix}\frac{d_{nm}-d_{nm}^0}{d_{nm}^0}\widehat{\textbf{r}}_{nm},&d_{nm}^0>0\\0,&d_{nm}^0\leq0\\\end{matrix}\right.&plus;\sum_{l,l\neq&space;m}\frac{d_{ml}^0-d_{ml}}{d_{ml}}\widehat{\textbf{r}}_{nm}\right)">,

where
* node <img src="https://latex.codecogs.com/svg.image?m"> is a neighbor of node <img src="https://latex.codecogs.com/svg.image?n">
* <img src="https://latex.codecogs.com/svg.image?\varsigma_{nm}"> is the stiffness of the bond between the nodes <img src="https://latex.codecogs.com/svg.image?n"> and <img src="https://latex.codecogs.com/svg.image?m"> (average of their individual stiffnesses)
* <img src="https://latex.codecogs.com/svg.image?d_{nm}"> is the distance between the nodes <img src="https://latex.codecogs.com/svg.image?n"> and <img src="https://latex.codecogs.com/svg.image?m">
* <img src="https://latex.codecogs.com/svg.image?d_{nm}^0"> is the equilibrium distance between the nodes <img src="https://latex.codecogs.com/svg.image?n"> and <img src="https://latex.codecogs.com/svg.image?m"> (average of their individual equilibrium distances)
* <img src="https://latex.codecogs.com/svg.image?\widehat{\textbf{r}}_{nm}"> is the unit vector in the direction of <img src="https://latex.codecogs.com/svg.image?\textbf{r}_m-\textbf{r}_n">
* node <img src="https://latex.codecogs.com/svg.image?l"> is another neighbor of node <img src="https://latex.codecogs.com/svg.image?n">
* <img src="https://latex.codecogs.com/svg.image?d_{ml}^0"> is the dimensionless equilibrium distance between the nodes <img src="https://latex.codecogs.com/svg.image?m"> and <img src="https://latex.codecogs.com/svg.image?l">; it is <img src="https://latex.codecogs.com/svg.image?2">, if node <img src="https://latex.codecogs.com/svg.image?l"> is opposite to node <img src="https://latex.codecogs.com/svg.image?m"> with respect to node <img src="https://latex.codecogs.com/svg.image?n">, and <img src="https://latex.codecogs.com/svg.image?\sqrt{2}"> otherwise
* <img src="https://latex.codecogs.com/svg.image?d_{ml}"> is the distance between the nodes <img src="https://latex.codecogs.com/svg.image?m"> and <img src="https://latex.codecogs.com/svg.image?l">

The first term inside the parentheses gives simple spring forces between node <img src="https://latex.codecogs.com/svg.image?n"> and its neighbors <img src="https://latex.codecogs.com/svg.image?m">, if both have greater-than-zero equilibrium length scales. If not, the term is omitted and the bond between the nodes can shrink and expand freely. The second term inside the parentheses (the sum over <img src="https://latex.codecogs.com/svg.image?l">) gives spring forces that try to keep the ratio of the distances between nodes <img src="https://latex.codecogs.com/svg.image?m"> and <img src="https://latex.codecogs.com/svg.image?l">, and <img src="https://latex.codecogs.com/svg.image?n"> and <img src="https://latex.codecogs.com/svg.image?m"> as 2 or <img src="https://latex.codecogs.com/svg.image?\sqrt{2}"> depending on if nodes <img src="https://latex.codecogs.com/svg.image?m"> and <img src="https://latex.codecogs.com/svg.image?l"> are opposite to each other with respect to node <img src="https://latex.codecogs.com/svg.image?n"> or not. Effectively this results in scale-independent angular forces that try to push the nodes into a square lattice. If a node is at the edge of the grid and does not have all four neighbors, the corresponding terms are simply omitted.

Note that there are no reaction forces; rather, rigid body motions and rotations can be eliminated in the relaxation. Alternatively, pixels can be fixed to prevent the system from drifting.

## Multiscale grid
To make relaxing the systems of nodes reasonably fast, the nodes are in a multiscale grid based on a full quadtree. The diagram below gives a one-dimensional example but this approach extends trivially to two dimensions.

<img src="https://github.com/petenez/deformer/blob/main/multiscale.png" width="800">

In this example, the red nodes have been fixed in place. The blue nodes are "active nodes" whose positions are updated during relaxation based on the forces they experience. An active node never has either active children or an active parent; rather, an active node must have either active nephews/nieces or active uncles/aunts, or both. The cyan and magenta nodes are "sampled nodes" whose positions are up- and downsampled from their parents below them and from their children above them, respectively. The sampled nodes simply provide couplings between the different levels of the multiscale grid. Neighbors for the active nodes to interact with are upsampled from their coarser parents' neighbors and downsampled from their neighbors' children. Downsampling simply takes the average of a node's children's positions. Upsampling on the other hand places a node a quarter of the way from its parent to its parent's neighbor. Given initial positions for the topmost "leaf nodes", all the other nodes' positions can be downsampled. After that, the active nodes positions can be updated and the sampled nodes up- and downsampled in an iterative fashion. When node positions are requested, the leaf nodes' positions can be upsampled starting from the active nodes. Note that when manipulating nodes (fixing or offsetting them, or setting their equilibrium length scale or their stiffness), the leaf nodes at the edges of the manipulated regions are treated. Inside and outside these regions, coarser grids of nodes are employed, as is demonstrated in the diagram above. For larger systems with a lower density of fixed and active leaf nodes, the reduction in computational effort can be orders of magnitude compared to a system where all nodes must be relaxed at full resolution.

Note also that I am familiar with multigrid methods but that is not what I wanted to implement here.

## Visualization
Deformed images are achieved by interpolating source image pixel positions from the leaf nodes' positions. The user needs to specify an origin and an offset vector in leaf node index coordinates. The origin indicates where the scanning for positions is started from and offset vectors define the step size from pixel to pixel in leaf node index coordinates. Index coordinates are simply the coordinate system set by the <img src="https://latex.codecogs.com/svg.image?i"> and <img src="https://latex.codecogs.com/svg.image?j"> indices of the nodes (but here we are not restricted to integer values). The offset vector for the second dimension of the image is calculated automatically (same magnitude and rotated by 90 degrees counterclockwise). The pixel positions are interpolated from the leaf nodes' positions using bicubic interpolation.

Finally, the target image can be formed from the source image pixels and their interpolated positions. The two methods of achieving this are either using a convolution or interpolating. The convolution method simply means convolving the pixels with a Gaussian kernel: Gaussian peaks are added to the target image according to a source image pixel's position. The peaks' amplitudes in the different color channels are proportional to the source image pixel's color. The user can choose the spread of the Gaussian kernel. The interpolation method is much more complicated because, instead of projecting each source image pixel onto the flat target image, one needs to solve the inverse problem of mapping each target image pixel to a position in the deformed source image. Having solved this problem, the rest is again simply interpolating the target image pixels from the source image using bicubic interpolation.

The main stages in the latter method are
* First iterate over the source image pixels. Find the target image pixel closest to each source image pixel's position and use Newton's method (starting from the source image pixel's indices) to find the position in source image index coordinates corresponding to the target image pixel's position. Interpolate the target image pixel's color using bicubic interpolation.
* Iterate over the target image pixels whose color is not set yet. Iterate over the neighbors of each unset pixel. If a neighbor is set, use Newton's method (starting from the position in source image pixel index coordinates corresponding to the set pixel's position) to find the position in source image index coordinates and interpolate the unset pixel's color. The unset pixel is now set. Keep iterating over the target image pixels in the same fashion while their number is decreasing.
* Possible remaining pixels' colors are obtained by diffusing from the neighboring set pixels.
