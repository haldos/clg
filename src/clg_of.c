/*
 * clg_of.c
 *
 * Copyright (C) 2012, Jorge Jara, DCC, SCIAN-Lab and BNI at Univ. of Chile
 *                     <jjara@dcc.uchile.cl>
 *                     Jose Delpiano, Univ. of the Andes <jdelpian@uandes.cl>
 *                     Steffen Haertel, SCIAN-Lab, BNI at Univ. of Chile
 *                     <shartel@med.uchile.cl>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "iio.c"
#include <math.h>   // Used for "exp"
#include <stdio.h>  // Used for "printf"
#include <stdlib.h> // Used for "free"

#define EPS 1E-12   // Precision threshold for linear system solving.
#define JROWS 3     // Size of array for storing J coeffs.
#define JCOLS 3
#define MIN_FILTER_SIZE_SIGMA 3
#define MIN_FILTER_SIZE_RHO   3

double **pmatrix(int nRows, int nCols);

void free_pmatrix(double **mat, int nRows, int nCols);

double lin2by2det(double a, double b, double c, double d);

void computeDerivatives(double **image1,
                        double **image2,
                        double **J[JROWS][JCOLS],
                        int nRows,
                        int nCols);

void matrixSmooth(double **matrix,
                  int nRows,
                  int nCols,
                  int filterSize,
                  double filterSigma);

void relax(double **u,
           double **v,
           double **J[JROWS][JCOLS],
           int nRows,
           int nCols,
           double alpha);

int calcCLG_OF(double* prevImage,
                double* currImage,
                double* uflow,
                double* vflow,
                int nCols,
                int nRows,
                int numIterations,
                double alpha,
                double rho,
                double sigma);

/**
 * pmatrix
 *
 * Allocates memory for a 2-dimensional array (matrix).
 * 
 * Parameters:
 * 
 * nRows The number of rows for the requested matrix (first dimension).
 * nCols The number of columns for the requested matrix (second dimension).
 */
double **pmatrix(int nRows, int nCols) {
	double **mat = (double **) malloc((size_t) (nRows * sizeof(double)));
	int i = 0;
    for (i=0; i<nRows; i++) {
        mat[i] = (double *) malloc((size_t) (nCols * sizeof(double)));
	}
	return mat;
}

/**
 * free_pmatrix
 * 
 * Free memory for a given 2-dimensional array (matrix).
 * 
 * Parameters:
 * 
 * mat   Pointer to the array to be freed.
 * nRows The number of rows for the requested matrix (first dimension).
 * nCols The number of columns for the requested matrix (second dimension).
 */
void free_pmatrix(double **mat, int nRows, int nCols) {
	int i = 0, j = 0;
    for (i=0; i<nRows; i++) {
		free(mat[i]);
	}
	free(mat);
}


/** 
 * lin2by2det
 *
 * Determinant for a 2x2 matrix, assuming the form
 * a b
 * c d
 */
inline double lin2by2det(double a, double b, double c, double d) {
    return a*d - b*c;
}


/**
 * computeDerivatives
 *
 * Discrete derivative computations for the input time frames (images). Spatial
 * derivatives are computed using the second image as the "current" time frame.
 *
 * Parameters:
 *
 * image1          Pointer to the first time frame (image).
 * image2          Pointer to the second time frame (image).
 * J[JROWS][JCOLS] Pointer to a motion tensor each frame (optical flow) pixel.
 *                 The motion tensor is a pointer to a matrix that will
 *                 contain the discrete derivatives. For a given pixel [i,j],
 *                 the tensor has the form
 *
 *                 J[0][0][i][j] = dfdx[i][j] * dfdx[i][j];
 *                 J[0][1][i][j] = dfdx[i][j] * dfdy[i][j];
 *                 J[0][2][i][j] = dfdx[i][j] * dfdt[i][j];
 *                 J[1][1][i][j] = dfdy[i][j] * dfdy[i][j];
 *                 J[1][2][i][j] = dfdy[i][j] * dfdt[i][j];
 *                 J[2][2][i][j] = dfdt[i][j] * dfdt[i][j];
 *
 *                 with f being the second (current) frame.
 *                 Since the tensor is symmetric, i.e.
 * 
 *                 J[1][0][i][j] = J[0][1][i][j];
 *                 J[2][0][i][j] = J[0][2][i][j];
 *                 J[2][1][i][j] = J[1][2][i][j];
 * 
 *                 only the upper half and the diagonal are stored.
 * nRows           Number of rows of the images.
 * nCols           Number of columns of the images.
 *
 */
void computeDerivatives(double **image1,
                        double **image2,
                        double **J[JROWS][JCOLS],
                        int nRows,
                        int nCols) {
    printf("  computeDerivatives\n");
    int i=1, j=1;

    // Image derivatives
    double **dfdx, **dfdy, **dfdt;
    dfdx = (double **) malloc((size_t) (nRows * sizeof(double)));
    dfdy = (double **) malloc((size_t) (nRows * sizeof(double)));
    dfdt = (double **) malloc((size_t) (nRows * sizeof(double)));

	// Using a 5x1 discrete derivative stencil according to Bruhn et al.
	// The operator derives from a 2nd order Taylor expansion.
    double h = 1.0 / nCols;
    double imDerStencil[] = {1.0 / 12 / h,
                            -8.0 / 12 / h,
                             0.0,
                             8.0 / 12 / h,
                            -1.0 / 12 / h};
    int imDerStencilSize = 5;
    int imDerStencilSide = imDerStencilSize / 2;

    printf("    loop for derivatives computation\n");
    for (i=0; i<nRows; i++) {

        dfdx[i] = (double *) malloc((size_t) (nCols * sizeof(double)));
        dfdy[i] = (double *) malloc((size_t) (nCols * sizeof(double)));
        dfdt[i] = (double *) malloc((size_t) (nCols * sizeof(double)));

        for (j=0; j<nCols; j++) {

            dfdx[i][j] = 0.0;
            dfdy[i][j] = 0.0;

            if (i>1 && i<nRows-2) {
                int count = 0;
                for (count=0; count<imDerStencilSize; count++) {
                    dfdx[i][j] += imDerStencil[count]
                        *image2[i+count-imDerStencilSide][j];
                }
            }
            if (j>1 && j<nCols-2) {
                int count = 0;
                for (count=0; count<imDerStencilSize; count++) {
                    dfdy[i][j] += imDerStencil[count]
                        *image2[i][j+count-imDerStencilSide];
                }
            }
            dfdt[i][j] = image2[i][j] - image1[i][j];
        }
    }

    i = 1;
    for (j=0; j<nCols; j++) {
        dfdx[i][j]         = (image2[2][j] - image2[0][j]) / 2.0 / h;
        dfdx[nRows-1-i][j] = (image2[nRows-1][j] - image2[nRows-3][j]) / 2.0 / h;
    }

    i = 0;
    for (j=0; j<nCols; j++) {
        dfdx[i][j] = dfdx[1][j];
        dfdx[nRows-i-1][j] = dfdx[nRows-2][j];
    }

    for (i=0;i<nRows;i++) {
        j = 1;
        dfdy[i][j] = (image2[i][2]-image2[i][0]) / 2.0 / h;
        dfdy[i][nCols-j-1] = (image2[i][nCols-1] - image2[i][nCols-3]) / 2.0 / h;
        j = 0;
        dfdy[i][j] = dfdy[i][1];
        dfdy[i][nCols-j-1] = dfdy[i][nCols-2];
    }

    printf("    storing J coeffs.\n");
    for (i=0; i<nRows; i++)
        for (j=0; j<nCols; j++) {
            J[0][0][i][j] = dfdx[i][j] * dfdx[i][j];
            J[0][1][i][j] = dfdx[i][j] * dfdy[i][j];
            J[0][2][i][j] = dfdx[i][j] * dfdt[i][j];
            J[1][1][i][j] = dfdy[i][j] * dfdy[i][j];
            J[1][2][i][j] = dfdy[i][j] * dfdt[i][j];
            J[2][2][i][j] = dfdt[i][j] * dfdt[i][j];
        }

    // Free memory
    for (i=0; i<nRows; i++) {
        free(dfdx[i]);
        free(dfdy[i]);
        free(dfdt[i]);
    }
    free(dfdx);
    free(dfdy);
    free(dfdt);
    printf("    done\n");
}


/**
 * relax
 *
 * Pointwise coupled Gauss-Seidel relaxation iteration for CLG-OF equations.
 * Each call to this function updates the current value of the solution,
 * u[1..m][1..n], v[1..m][1..n], using the motion tensor
 * J[JROWS][JCOLS][1..m][1..n].
 * Neumann boundary conditions are used (derivatives are set to zero).
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly 
 *                 smoothed) derivatives for each image/optical flow pixel.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 * alpha           Optical flow global smoothing coefficient.
 */
void relax(double **u,
           double **v,
           double **J[JROWS][JCOLS],
           int nRows,
           int nCols,
           double alpha) {

    int i, ipass, isw, j, jsw=1;
    double hx, hy, hx2, hy2, discr, detU, detV;

    hx = 1.0 / (nCols-1);
    hy = hx;
    hx2 = hx * hx;
    hy2 = hy * hy;

    for (ipass=1; ipass<=2; ipass++, jsw=3-jsw) {       //Red and black sweeps.
        isw = jsw;
        for (j=0; j<=nCols-1; j++, isw=3-isw)
            for (i=isw-1; i<=nRows-1; i+=2) {           //Gauss-Seidel formula.
                if (i>0 && i<nRows-1) {
                    if (j>0 && j<nCols-1) {
                        discr = lin2by2det(
                                 alpha/hx2*2.0 + alpha/hy2*2.0 + J[0][0][i][j],
                                 J[0][1][i][j],
                                 J[0][1][i][j],
                                 alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j]);

                        if (abs((double)discr) > EPS) {
                            detU = lin2by2det(
                                     alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                       + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                       - J[0][2][i][j],
                                     J[0][1][i][j],
                                     alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                       + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                       - J[1][2][i][j],
                                     alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                       + J[1][1][i][j]);

                            u[i][j]= detU / discr;

                            detV = lin2by2det(
                                     alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                       + J[0][0][i][j],
                                     alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                       + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                       - J[0][2][i][j],
                                     J[0][1][i][j],
                                     alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                       + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                       - J[1][2][i][j]);

                            v[i][j] = detV / discr;

                        } else {
                            u[i][j] =  (alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                      + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                      - (J[0][1][i][j]*v[i][j] + J[0][2][i][j]))
                                    / (alpha / hx2*2.0 + alpha / hy2*2.0
                                      + J[0][0][i][j]);

                            v[i][j] =  (alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                      + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                      - (J[0][1][i][j]*u[i][j] + J[1][2][i][j]))
                                    / (alpha / hx2*2.0 + alpha / hy2*2.0
                                      + J[1][1][i][j]);
                        }
                    } else if (j==0) {
                        u[i][j] = u[i][j+1];
                        v[i][j] = v[i][j+1];
                    } else if (j==nCols-1) {
                        u[i][j] = u[i][j-1];
                        v[i][j] = v[i][j-1];
                    } else {
                        printf("  index error(2).\n");
                        return;
                    }
                } else if (i==0) {
                    if (j>0 && j<nCols-1) {
                        u[i][j] = u[i+1][j];
                        v[i][j] = v[i+1][j];
                    } else if (j==0) {
                        u[i][j] = u[i+1][j+1];
                        v[i][j] = v[i+1][j+1];
                    } else if (j==nCols-1) {
                        u[i][j] = u[i+1][j-1];
                        v[i][j] = v[i+1][j-1];
                    } else {
                        printf("  index error.\n");
                        return;
                    }
                } else if (i==nRows-1) {
                    if (j>0 && j<nCols-1) {
                        u[i][j] = u[i-1][j];
                        v[i][j] = v[i-1][j];
                    } else if (j==0) {
                        u[i][j] = u[i-1][j+1];
                        v[i][j] = v[i-1][j+1];
                    } else if (j==nCols-1) {
                        u[i][j] = u[i-1][j-1];
                        v[i][j] = v[i-1][j-1];
                    } else {
                        printf("  index error(3).\n");
                        return;
                    }
                }
            }
    }
}


/**
 * matrixSmooth
 *
 * Performs a gaussian smoothing on a given matrix, with a specified size for
 * the kernel (filter) and standard deviation value.
 *
 * Parameters:
 * nrMatrix     Pointer to the input matrix.
 * nRows        Number of rows of the input matrix.
 * nCols        Number of columns of the input matrix.
 * filterSize   Size of the gaussian kernel to be computed and applied.
 * filterSigma  Standard deviation for the gaussian kernel.
 */
void matrixSmooth(double **nrMatrix,
                  int rows,
                  int cols,
                  int filterSize,
                  double filterSigma) {

    int i=1, j=1, k=0, l=0;
    double sum = 0.0;
    double **aux;

    printf("  matrixSmooth\n");
    printf("    allocating aux. matrix\n");
    // Matrix backup
    aux = (double **) malloc((size_t) rows * sizeof(double));
    for (i=0; i<rows; i++) {
        aux[i] = (double *) malloc((size_t) cols * sizeof(double));
        for (j=0; j<cols; j++) {
            aux[i][j] = nrMatrix[i][j];
        }
    }

    // Gaussian filter 1D
    printf("    creating gaussian filter 1D\n");

    int filterSide = filterSize / 2; // int. division, odd filterSize assumed
    double *filter = (double *) malloc((size_t) filterSize * sizeof(double));

    for (k=0; k<filterSize; k++) {
        //filter = exp(-(x-N-1).^2/2/sigma^2);
        filter[k] = exp(-(k-filterSide) * (k-filterSide)
                        / 2.0 / filterSigma / filterSigma);
        sum += filter[k];
    }

    // Convolution
    double momentSum = 0.0;
    sum = 0.0;
    int currentPos = 0;

    printf("    convolution, 1st dim\n");

    // First dimension
    
    // non-boundary pixels
    for (i=filterSide; i<rows-filterSide; i++) {
        for (j=0; j<cols; j++) {
            nrMatrix[i][j] = 0.0;
            for (k=0; k<filterSize; k++) {
                nrMatrix[i][j] += aux[i+k-filterSide][j] * filter[k];
            }
        }
	}

    // upper boundary
    for (i=0; i<filterSide; i++) {
        for (j=1; j<cols; j++) {
            sum = 0.0;
            momentSum = 0.0;
            for (k=0; k<filterSize; k++) {
                currentPos = i+k-filterSide;
                if (currentPos>=0 && currentPos<rows) {
                    momentSum += aux[currentPos][j] * filter[k];
                    sum += filter[k];
                }
            }
            nrMatrix[i][j] = momentSum / sum;
        }
    }

    // lower boundary
    for (i=rows-1-filterSide; i<rows; i++) {
        for (j=0; j<cols; j++) {
            sum = 0.0;
            momentSum = 0.0;
            for (k=0; k<filterSize; k++) {
                currentPos = i+k-filterSide;
                if (currentPos>=0 && currentPos<rows) {
                    momentSum += aux[currentPos][j] * filter[k];
                    sum += filter[k];
                }
            }
            nrMatrix[i][j] = momentSum / sum;
        }
    }

    // Matrix backup
    printf("    1st dim, matrix backup\n");
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            aux[i][j] = nrMatrix[i][j];
        }
    }

    printf("    convolution, 2nd dim\n");
    // Second dimension

	// central pixels
    for (i=0; i<rows; i++) {
        for (j=filterSide; j<cols-filterSide; j++) {
            nrMatrix[i][j] = 0.0;
            for (k=0; k<filterSize; k++) {
                nrMatrix[i][j] += aux[i][j+k-filterSide] * filter[k];
            }
        }
	}

	// left boundary
    for (i=0; i<rows; i++) {
        for (j=0; j<filterSide; j++) {
            sum = 0.0;
            momentSum = 0.0;
            for (k=0; k<filterSize; k++) {
                currentPos = j+k-filterSide;
                if (currentPos>=0 && currentPos<cols) {
                    momentSum += aux[i][currentPos] * filter[k];
                    sum += filter[k];
                }
            }
            nrMatrix[i][j] = momentSum / sum;
        }
    }

	// right boundary
    for (i=0; i<rows; i++) {
        for (j=cols-1-filterSide; j<cols; j++) {
            sum = 0.0;
            momentSum = 0.0;
            for (k=0; k<filterSize; k++) {
                currentPos = j+k-filterSide;
                if (currentPos>=0 && currentPos<cols) {
                    momentSum += aux[i][currentPos] * filter[k];
                    sum += filter[k];
                }
            }
            nrMatrix[i][j] = momentSum / sum;
        }
    }

    // Free memory
    for (i=0; i<rows; i++)
        free(aux[i]);
    free(filter);
    free(aux);
    printf("    done\n");
}


/**
 * calcCLG_OF
 *
 * Main CLG-optical flow (CLG-OF) computation function.
 *
 * Parameters:
 *
 * prevImage     Pointer to the first image (the "previous" time frame).
 * currImage     Pointer to the second image (the "current" time frame).
 * uflow         Pointer to the horizontal component of the CLG-OF solution.
 * vflow         Pointer to the vertical component of the CLG-OF solution.
 * nCols         Number of image columns (same for the CLG-OF vector field).
 * nRows         Number of image rows (same for the CLG-OF vector field).
 * numIterations Number of iterations for iterative solution.
 * alpha         Global smoothing coefficient of the CLG-OF.
 * rho           Local spatio-temporal smoothing coefficient of the CLG-OF.
 * sigma         Standard deviation for the optional gaussian smoothing kernel,
 *               applied to the input images prior to the CLG-OF computation.
 *
 */
int calcCLG_OF(double* prevImage,
               double* currImage,
               double* uflow,
               double* vflow,
               int nCols,
               int nRows,
               int numIterations,
               double alpha,
               double rho,
               double sigma) {

    printf("calc_clg\n");
    printf("  setting up variables\n");
    int i=0, j=0;

    // IDL variables to NR variables
    double **prevFrame, **currFrame;
    prevFrame = pmatrix(nRows, nCols);
    currFrame = pmatrix(nRows, nCols);

    for (j=0; j<nRows; j++) {
        for (i=0; i<nCols; i++) {
            prevFrame[j][i] = prevImage[i + j*nCols];
            currFrame[j][i] = currImage[i + j*nCols];
        }
    }

    printf("  allocating memory for arrays\n");
    double **u, **v;
    u = pmatrix(nRows, nCols);
    v = pmatrix(nRows, nCols);

    for (j=0; j<nRows; j++) {
        for (i=0; i<nCols; i++) {
            u[j][i] = uflow[i + j*nCols];
            v[j][i] = vflow[i + j*nCols];
        }
    }

    printf("  allocating memory for derivatives matrices\n");
    double **J[JROWS][JCOLS];
    int k=0, l=0;
    for (k=0; k<JROWS; k++)
        for (l=0; l<JCOLS; l++)
            J[k][l] = pmatrix(nRows, nCols);

    i=1, j=1;

    // Gaussian pre smoothing
    int filterSize = (int) (2.0 * (int)(2.5 * sigma) + 1.0);
    printf("  filter size %i for sigma %f\n", filterSize, sigma);

    if (filterSize >= MIN_FILTER_SIZE_SIGMA) {
        printf("  apply gaussian smoothing for each frame\n");
        matrixSmooth(prevFrame, nRows, nCols, filterSize, sigma);
        matrixSmooth(currFrame, nRows, nCols, filterSize, sigma);
    } else {
		printf("  no gaussian smoothing was applied to input frames\n");
	}

    // J coeffs obtention
    computeDerivatives(prevFrame, currFrame, J, nRows, nCols);

    // Odd filter size
    filterSize = 2 * ((int) (2.5 * rho)) + 1;
    printf("  filter size %i for rho %f\n", filterSize, rho);

    if (filterSize < MIN_FILTER_SIZE_RHO)
        filterSize = MIN_FILTER_SIZE_RHO;

    k=0, l=0;
    if (filterSize >= 2) {
        printf("  applying local spatio temporal smoothing\n");
        for (k=0; k<JROWS; k++)
            for (l=0; l<JCOLS; l++)
                matrixSmooth(J[k][l], nRows, nCols, filterSize, rho);
    } else {
        printf("  no local spatio temporal smoothing was applied\n");
	}

    if (numIterations == 0)
        numIterations = (int)(nRows * nCols / 8.0);

    printf("  performing %i relax iterations\n", numIterations);
    int count = 0;
    for (count=0; count<numIterations; count++)
			relax(u, v, J, nRows, nCols, alpha);

    // Optical flow scaling (image supposed to be rx1.0)
    printf("  scaling output vectors u,v \n");
    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            u[i][j] = nCols * u[i][j];
            v[i][j] = nCols * v[i][j];
        }
    }

    printf("  formatting output\n");
    for (j=0; j<nRows; j++) {
        for (i=0; i<nCols; i++) {
            uflow[i + j*nCols] = u[j][i];
            vflow[i + j*nCols] = v[j][i];
        }
    }

    // Free memory
    printf("  freeing memory\n");
    free_pmatrix(u, nRows, nCols);
    free_pmatrix(v, nRows, nCols);

    for (k=0; k<JROWS; k++) {
        for (l=0; l<JCOLS; l++) {
            free_pmatrix(J[k][l], nRows, nCols);
        }
    }
    free_pmatrix(prevFrame, nRows, nCols);
    free_pmatrix(currFrame, nRows, nCols);

    printf("calc_clg: done\n");
    return 1;
}


/**
 * Demo for CLG OF computation, main function.
 * Image loading
 * Call of function calc_clg
 * Flow file saving
 */
int main(int argc, char *argv[]) {
    if (argc != 8) {
        printf("Usage: %s input_image_1 input_image_2 alpha rho sigma numIt output\n", argv[0]);
    } else {
        int w, h, pixeldim;

        // Load specified images
        float *im1 = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
        float *im2 = iio_read_image_float_vec(argv[2], &w, &h, &pixeldim);
        printf("Two images loaded:\n\tim,: %dx%d image with %d channel(s)\n\tim2, %dx%d image with %d channel(s)\n", w, h, pixeldim, w, h, pixeldim);

		// read parameters from console input
		double alpha = (double) atof(argv[3]);
		double rho   = (double) atof(argv[4]);
		double sigma = (double) atof(argv[5]);
		int numIt    = atoi(argv[6]);

        // CLG optical flow computation
        printf("allocating memory for data structures\n");
        double* u;
        double* v;
        double* i1;
        double* i2;
        u = (double *) malloc((size_t) h * w * sizeof(double));
        v = (double *) malloc((size_t) h * w * sizeof(double));
        i1 = (double *) malloc((size_t) h * w * sizeof(double));
        i2 = (double *) malloc((size_t) h * w * sizeof(double));

        printf("copying and casting data for double precision\n");
        int i = 0;
        int j = 0;
        for (i=0; i<w; i++) {
            for (j=0; j<h; j++) {
                u[w*j+i] = 0.0;
                v[w*j+i] = 0.0;
                i1[w*j+i] = ((double) im1[w*j+i]);
                i2[w*j+i] = ((double) im2[w*j+i]);
            }
        }

        printf("call to calc_clg\n");
        int res = calcCLG_OF(i2, i1, u, v, w, h, numIt, alpha, rho, sigma);

        // Save image

        float* uf;
        float* vf;
        uf = (float *) malloc((size_t) h * w * sizeof(float));
        vf = (float *) malloc((size_t) h * w * sizeof(float));

        float *x = (float*)malloc(w * h * 2 * sizeof * x);
        for (i=0; i<w; i++) {
            for (j=0; j<h; j++) {
                x[2*(w*j+i)]   = u[w*j+i];
                x[2*(w*j+i)+1] = v[w*j+i];

                uf[w*j+i] = u[w*j+i];
                vf[w*j+i] = v[w*j+i];
            }
        }

        iio_save_image_float_vec(argv[7], x, w, h, 2);
        iio_save_image_float_vec("flow_u.png", uf, w, h, 1);
        iio_save_image_float_vec("flow_v.png", vf, w, h, 1);

        //free memory
        free(im1);
        free(im2);
        free(x);
        free(uf);
        free(vf);
        free(u);
        free(v);
        printf("clg-optical flow computation done.\n");
        return 0;
    }
}
