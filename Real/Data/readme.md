The original dataset Xsus.mat is too large to be uploaded directly. It is in the shape of 3001 * 100 * 100, where the 1st dimension is Fourier coefficients, and the 2-3 dimensions are for spatial axes.

The Xsus.mat is divided into 10 subfiles. Each file is 300 * 100 * 100 size, except for the final subfile 301 * 100 * 100. They can be concatenated along the 1st dimension according to the indices to rebuild the full Xsus.mat.
