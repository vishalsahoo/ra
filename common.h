#include <vector>
#include <cstdlib>
#include <bitset>
#include <map>            //std::map
#include <cstring>        // memcpy, memset
//#include "bmmatrix.cpp"
#include "bitset.h"
// --- Problem parameters ---
// NX can be a power of 2, min val 8
#define NX 8

#define NY NX
#define DIM 2
// For a diamond shape with index ranges from (1,NX-2)x(1,NY-2)
#define NEB ((NX-2)+(NY-2))
#define NTOT (NX*NY)
#define BLOCKSIZE 64
#define NR 8
#define NC 8
// Covers the whole domain

struct arr_t
{
  double data[NY][NX]; // NB: [j][i] indexing
  double* ptr=&(data[0][0]); // NB: assume allocated in contiguous block!
  // NB: For large NX this might break
};

// Struct for holding eb index data
struct ebix_t
{
  int ebix[NEB][DIM]; // ij indices of EB cells
  int inout[NY][NX]; // flag over entire domain, out=-1, eb=0, irreg=1, reg=2
};

// Struct for holding eb sparse matrix information
struct ebsparseop_t
{
  std::vector<size_t> index; // non-zero global offset index for result
  std::vector<std::vector<size_t> > offset; // non-zero global offset indices
  std::vector<std::vector<double> > weight; // non-zero values for each index
  std::vector<double> bc; // values to add for bc's
};
// created x and y vectors for bcsr computation in cuda
struct ebvector{
    double *data;
    int size = 0;
};

struct ebsparsematrix_t
{
  size_t N=NTOT; // size of matrix, N x N
  std::vector<size_t> nzrow; // for each non-zero, the row / global index
  std::vector<size_t> nzcol; // for each non-zero, the col index
  std::vector<double> entry; // non-zero values for each index
};

struct block_t {
    double matrix[BLOCKSIZE][BLOCKSIZE] = {{ 0 }}; // A dense 0 padded matrix of the non-zero values
    size_t row; // The starting row of the block
    size_t col; // The starting col of the block
};

struct ebbcsrmatrix_t
{
    size_t blocksize = BLOCKSIZE; // Size of the blocks B*B
    std::vector<double> values; // The vector of blocks
    std::vector<size_t> cols; //vector of columns for nnz blocks. The size of this is the num of blocks
    std::vector<size_t> block_row_ptr; // Compressed block row pointer
};

struct ebcsrmatrix_t {
    std::vector<double> values; // The vector of values
    std::vector<size_t> cols; //vector of columns for nnz blocks
    std::vector<size_t> row_ptr; // Compressed row pointer
};

struct ebbcsrBMmatrix_t
{
std::vector<double> values;
std::vector<size_t> cols;
std::vector<size_t> block_row_ptr;
std::vector<size_t> value_pos;
std::vector<BitSet<BLOCKSIZE, BLOCKSIZE> > bit;
};

/**
 * @brief Converting from COO format to BCSR format
 * Jake + Priya
 * @param ebmat COO
 * @param ebbcsr BCSR
 */
void convertToBCSR(ebsparsematrix_t& ebmat, ebbcsrmatrix_t& ebbcsr);

/**
 * @brief Converting from COO format to CSR format.
 * Jake
 * @param ebmat COO
 * @param ebcsr CSR
 */
void cootocsr(ebsparsematrix_t& ebmat, ebcsrmatrix_t& ebcsr);

// Solve the bcsr for the given rhs (which should include bc's)
// NB: note, maxiter should O(NX^2) to get to convergence
void solvebcsr(const ebbcsrmatrix_t& ebbcsr, const ebsparseop_t& ebop,
    arr_t& soln, arr_t& rhs, arr_t& tmp, int maxiter, double tol);
//
//void applybcsr(const ebbcsrmatrix_t& ebbcsr, const arr_t& soln, arr_t& opval);
//

/**
 * @brief Calling global cuda kernel to perform computation on bcsr.
 * Priya
 * @param ebbcsr bcsr
 * @param x multiplier vector
 * returns a vector that includes one element per row after performing x*bcsr computation
 */
ebvector spmvbcsrcu(ebbcsrmatrix_t& ebbcsr, ebvector& x);

/**
 * @brief Calculate time required to perform calculations in sparse and dense block
 * @param v A BCSR block
 * @param sparse_compact_row Vector containing values without zero padding
 * @param bitmap Bitmap of blocks with values
 * @param block_row_ptr Block row data
 */

std::vector<double> calculateTime(ebbcsrBMmatrix_t bmp, std::vector<double> multiplier);

/**
 * @brief Generate bitmap and sparse compact row i.e. values without zero padding
 * @param ebbcsrmatrix structure
 * @param ebbcsrBMmatrix_t structure
 */
//void  generateBitmap(std::vector<block_t>& v, std::vector<double> &sparse_compact_row, std::vector<std::bitset<BLOCKSIZE*BLOCKSIZE>>& bitmap);
void generateBitmap(ebbcsrmatrix_t& ebbcsr, ebbcsrBMmatrix_t &bmp);

/**
 * @brief Decompress BCSR into original matrix with zero padding
 * @param bitmap structure
 */
//void decompress(std::vector<double> sparse_compact_row, std::vector<std::bitset<BLOCKSIZE*BLOCKSIZE>> bitmap, std::vector<size_t> block_row_ptr, std::vector<size_t> columns);
void decompress(ebbcsrBMmatrix_t bmp, ebsparsematrix_t eb);

/**
 * @brief Call cuda kernel
 * @param bitmap structure
 * @param x multiplier vector
 */
std::vector<double> compactvaluescu(ebbcsrBMmatrix_t& bmp, std::vector<double>& x);


