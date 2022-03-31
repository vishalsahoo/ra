// ----------------------------------------------------------------------------
// Test problem for demonstrating a sparse embedded boundary (EB)
// Laplacian operator in a diamond geometry on a single domain 
// ----------------------------------------------------------------------------

#include <iostream>       // std::cout
#include <string>         // std::string
#include <bitset>         // std::bitset
#include <tuple>          // std::tuple
#include <vector>         // std::vector
#include <cstring>        // memcpy, memset
#include <cassert>        // assert
#include <map>            // std::map
#include <chrono>         // std::chrono
//#include "common.h"
#include "timelog.cpp"
#include "bitset.h"
#include "coo2bcsr.cpp" //for g++ compiler TO BE DELETED FOR MAKE
#include "bmmatrix.cpp" //for g++ compiler TO BE DELETED FOR MAKE

// Sizes for arrays and loops
const size_t dim[DIM] = {NX, NY};
const int lo[DIM] = {0, 0}; // min index
const int hi[DIM] = {NX-1, NY-1}; // max index
const size_t neb = NEB; // number of eb partial cells
const double dx = 1./(double) NX; // grid spacing
const int bsize = BLOCKSIZE;

// all-index memset to =0
void zero(arr_t& arr)
{
  memset(arr.ptr, 0, NX*NX*sizeof(double));
}

// all-index memcpy: from->to
void copy(arr_t& to, const arr_t& from)
{
  memcpy(to.ptr, from.ptr, NX*NX*sizeof(double));
}

// Our exact solution = x^2 + y^2, so Laplacian = 4
double exact(const int& i, const int& j)
{
  double x = (.5 + i)*dx;
  double y = (.5 + j)*dx;
  return (x*x+y*y);
}

// regular cell, 5pt Laplacian standard stencil
void stencil_reg(int i, int j, std::vector<size_t>& offset, 
    std::vector<double>& weight)
{
  //     1
  // 1 *-4*  1   * means center
  //     1
  offset.push_back(i + dim[0]*j); // index offset of this cell
  weight.push_back(-4); // weight of this cell
  for (int s=0; s < 4; ++s)
  {
    int di = (s/2)*(2*(s%2)-1); // i offset, -1,1,0,0
    int dj = (1-s/2)*(2*(s%2)-1); // j offset, 0,0,-1,1
    offset.push_back(i+di + dim[0]*(j+dj)); // valid index
    weight.push_back(1); // weight=1, standard 5pt
  }
}

// eb stencil, using only inside points with bc corrections
void stencil_eb(int i, int j, ebix_t& eb, std::vector<size_t>& offset, 
    std::vector<double>& weight, double& bc)
{
  //      1
  // bc *-4*  1   * means center
  //     bc 
  
  // Don't need this test
  // bool corner = (abs(i-NX/2+.5)==.5 || abs(j-NY/2+.5)==.5);
  offset.push_back(i + dim[0]*j); // index offset of this cell
  weight.push_back(-4); // weight of this cell
  bc = 0;
  // for each of the other points, if outside add the bc value
  for (int s=0; s < 4; ++s)
  {
    int di = (s/2)*(2*(s%2)-1); // i offset, -1,1,0,0
    int dj = (1-s/2)*(2*(s%2)-1); // j offset, 0,0,-1,1
    if (eb.inout[j+dj][i+di] < 0)
    {
      bc += exact(i+di,j+dj);
    }
    else
    {
      offset.push_back(i+di + dim[0]*(j+dj)); // valid index
      weight.push_back(1); // weight=1, standard 5pt
    }
  }
}

// Initialize the solution, eb index data, and eb sparse matrix
void initialize(arr_t& soln, ebix_t& eb, ebsparseop_t& ebop, 
    ebsparsematrix_t& ebmat)
{
  // Initialize the mask of eb in/out data for the diamond
  // std::cout << "Setup for in/out: " << std::endl;
  size_t n=0;
  for (int j=lo[1]; j <= hi[1]; ++j)
  for (int i=lo[0]; i <= hi[0]; ++i)
  {
    const float di = abs(NX - 2*i - 1);
    const float dj = abs(NY - 2*j - 1);
    int inout = (di+dj <= NX-2) - (di+dj >= NX-2); 
    // Test to determine if [j][i] is on the EB
    // bool oneb = ((i+j)==(NX/2)) || ((i-j)==(NX/2-1))
    //             || ((j-i)==(NX/2-1)) || ((i+j)==(3*NX/2-2));
    // oneb = oneb && (i>lo[0] && i<hi[0] && j>lo[1] && j<hi[0]);
    // std::cout << "  [" << j << "][" << i << "]= " 
    //   << ", in/out/val = " << inout << std::endl;
    // ", on eb = " << oneb << std::endl;
    eb.inout[j][i] = inout;
    if (inout==0)
    {
      eb.ebix[n][0] = i;
      eb.ebix[n][1] = j;
      // std::cout << "EB[" << n << "] = [" << j << "][" << i << "]= " 
      //  << ", in/out = " << inout << std::endl;
      ++n;
    }
  }
  assert(n == neb);

  // Initialize the soln, op to cell center values
  // soln = x^2 + y^2
  // op = (d_xx + d_yy) soln = 4
  for (int j=lo[1]; j <= hi[1]; ++j)
  for (int i=lo[0]; i <= hi[0]; ++i)
  {
    if (eb.inout[j][i] < 0) // outside domain
      soln.data[j][i] = 0;
    else
      soln.data[j][i] = exact(i,j);
  }

  // Initialize the eb sparse operator and matrix
  int nnz=0;
  for (int j=lo[1]; j <= hi[1]; ++j)
  for (int i=lo[0]; i <= hi[0]; ++i)
  {
    const size_t index = i+dim[0]*j;
    std::vector<size_t> offset;
    std::vector<double> weight;
    double bc=0;
    if (eb.inout[j][i]==-1)
      continue;
    if (eb.inout[j][i]>0)
      stencil_reg(i, j, offset, weight);
    else
      stencil_eb(i, j, eb, offset, weight, bc);

    ebop.index.push_back(index);
    ebop.offset.push_back(offset);
    ebop.weight.push_back(weight);
    ebop.bc.push_back(bc);

    for (int k=0; k < offset.size(); ++k)
    {
      ebmat.nzrow.push_back(index);
      ebmat.nzcol.push_back(offset[k]);
      ebmat.entry.push_back(weight[k]);
      std::cout<<std::endl<<index<<","<<offset[k]<<","<<weight[k];
      ++nnz;
    }
  }
  assert(nnz = ebmat.entry.size());
}


// Prints out the contents of the array for [0,NX)x [0,NY)
void dump(arr_t& arr, std::string msg)
{
  std::cout << "Values of " << msg << " in cell [j][i]:" << std::endl;
  for (int j=0; j < dim[1]; ++j)
  for (int i=0; i < dim[0]; ++i)
  {
    std::cout << "  [" << j << "][" << i << "]= " 
      << arr.data[j][i] << std::endl;
  }
}

// Prints out the contents of the eb sparse operator
void dumpop(ebix_t& eb, ebsparseop_t& ebop)
{
  std::cout << "Sparse op offsets, values:" << std::endl;
  for (int n=0; n < ebop.index.size(); ++n)
  {
    int i = ebop.index[n] % dim[0];
    int j = (ebop.index[n] - i)/dim[0];
    const std::vector<size_t>& offset = ebop.offset[n];
    const std::vector<double>& weight = ebop.weight[n];

    std::cout << "  For index [" << j << "][" << i << "]: " << std::endl
      << "    bc value: " << ebop.bc[n] << std::endl;
    for (int ix=0; ix < offset.size(); ++ix)
    {
      int di = offset[ix] % dim[0] - i;
      int dj = (offset[ix] - di)/dim[0] - j;
      std::cout << "    [j" << ((dj<0)?"":"+") << dj 
        << "][i" << ((di<0)?"":"+") << di << "]: " << weight[ix] << std::endl;
    }
  }
}

// Prints out the contents of the eb sparse matrix
// std=true means a std::cout format, while false means matlab
void dumpmat(ebsparsematrix_t& ebmat, bool std=true)
{
  size_t N = ebmat.N;
  if (std)
  {
    std::cout << "Sparse matrix size: " << N << " x " << N << std::endl;
    std::cout << "  Nonzeros, [row,col]=entry" << std::endl;
  }
  else // matlab
    std::cout << "A = sparse(" << N << "," << N << ");" << std::endl;
  for (int n=0; n < ebmat.entry.size(); ++n)
  {
    int i = ebmat.nzrow[n];
    int j = ebmat.nzcol[n];
    double e = ebmat.entry[n];
    if (std)
      std::cout << "  [" << i << "][" << j << "] = " << e << std::endl;
    else // matlab
    std::cout << " A(" << i+1 << "," << j+1 << ") = " << e << ";" << std::endl;
  }
}

// Apply the eb sparse operator to soln, store result in opval
void applyop(const ebsparseop_t& ebop, const arr_t& soln, 
    arr_t& opval, bool use_bc)
{
  const int ncell = ebop.index.size();
  for (int n=0; n < ncell; ++n)
  {
    const size_t& ix = ebop.index[n];
    *(opval.ptr+ix) = (use_bc) ? ebop.bc[n] : 0;
    const std::vector<size_t>& offset = ebop.offset[n];
    const std::vector<double>& weight = ebop.weight[n];
    for (int s=0; s < offset.size(); ++s)
      *(opval.ptr+ix) += *(soln.ptr+offset[s])*weight[s];
  }
}

// Apply the eb sparse matrix to soln, store result in opval
void applymat(const ebsparsematrix_t& ebmat, const arr_t& soln, arr_t& opval)
{
  size_t N = ebmat.N;
  for (int n=0; n < ebmat.entry.size(); ++n)
  {
    const int i = ebmat.nzrow[n];
    const int j = ebmat.nzcol[n];
    const double e = ebmat.entry[n];
    opval.ptr[i] += soln.ptr[j]*e;
  }
}

// Apply the eb BCSR to soln, store result in opval
void applybcsr(const ebbcsrmatrix_t& ebbcsr, const arr_t& soln, arr_t& opval) {
  int blocksize = ebbcsr.blocksize; // Small optimization to move this into a local variable

  for (int i = 0; i < ebbcsr.block_row_ptr.size() - 1; i++) {
    int rowstart = ebbcsr.block_row_ptr[i];
    int rowend = ebbcsr.block_row_ptr[i + 1];

    // This loop goes over every block in a non-zero row
    for (int block = rowstart; block < rowend; block++) {
      int column = ebbcsr.cols[block] * blocksize; // Significant optimization by keeping this multiply out of the a += x*y operation
      // These loops go through the block
      for (int ii = 0; ii < BLOCKSIZE; ii++) {
        for (int jj = 0; jj < BLOCKSIZE; jj++) {
          opval.ptr[blocksize * i + ii] += soln.ptr[column + jj] * ebbcsr.values[block * blocksize * blocksize + ii * blocksize + jj];
        }
      }
    }
  }
}

// Apply the boundary condition correction to the vector
void applybc(const ebsparseop_t& ebop, arr_t& tmp, const double scale)
{
  const int ncell = ebop.index.size();
  for (int n=0; n < ncell; ++n)
  {
    const size_t& ix = ebop.index[n];
    tmp.ptr[ix] += scale*ebop.bc[n];
  }
}

// Solve the eb sparse operator w/ bc's for the given rhs
// NB: note, maxiter should O(NX^2) to get to convergence
void solve(const ebsparseop_t& ebop, arr_t& soln, arr_t& rhs, 
    arr_t& tmp, int maxiter, double tol)
{
  std::cout << "Solver: maxiter = " << maxiter << ", tolerance = " 
    << tol << std::endl;
  const int ncell = ebop.index.size();
  int iter=0;
  double resid=tol;
  double param=.125; // =1/8, only valid for this problem
  // Correct the rhs for the bc values
  for (int n=0; n < ncell; ++n)
    *(rhs.ptr+ebop.index[n]) -= ebop.bc[n];

  std::cout << "   beginning solver point Jacobi iterations ..." << std::endl;
  double maxresid;
  while (iter < maxiter)
  {
    applyop(ebop, soln, tmp, false);
    // Point Jacobi update
    maxresid=0;
    for (int n=0; n < ncell; ++n)
    {
      const size_t ix=ebop.index[n];
      const double resid = *(tmp.ptr+ix) - *(rhs.ptr+ix);
      maxresid = std::max<double>(maxresid, abs(resid));
      *(soln.ptr+ix) += resid*param;
    }

    if (maxresid < tol)
      break;
    ++iter;
  }
  std::cout << "Solver exited, " << iter << " iterations, residual = "
    << maxresid << std::endl;
}

// Solve the eb sparse matrix for the given rhs (which should include bc's)
// NB: note, maxiter should O(NX^2) to get to convergence
void solvesparse(const ebsparsematrix_t& ebmat, const ebsparseop_t& ebop, 
    arr_t& soln, arr_t& rhs, arr_t& tmp, int maxiter, double tol)
{
  //std::cout << "Solver: maxiter = " << maxiter << ", tolerance = " 
  //  << tol << std::endl;
  const int ncell = ebop.index.size();
  int iter=0;
  double resid=tol;
  double param=.125; // =1/8, only valid for this problem

  //std::cout << "   beginning solver point Jacobi iterations ..." << std::endl;
  double maxresid;
  while (iter < maxiter)
  {
    zero(tmp);
    applymat(ebmat, soln, tmp);
    // Point Jacobi update
    maxresid=0;
    for (int n=0; n < ncell; ++n)
    {
      const size_t ix=ebop.index[n];
      const double resid = tmp.ptr[ix] - rhs.ptr[ix];
      maxresid = std::max<double>(maxresid, abs(resid));
      soln.ptr[ix] += resid*param;
    }

    if (maxresid < tol)
      break;
    ++iter;
  }
  //std::cout << "Solver exited, " << iter << " iterations, residual = "
  //  << maxresid << std::endl;
}

// Solve the bcsr for the given rhs (which should include bc's)
// NB: note, maxiter should O(NX^2) to get to convergence
void solvebcsr(const ebbcsrmatrix_t& ebbcsr, const ebsparseop_t& ebop,
    arr_t& soln, arr_t& rhs, arr_t& tmp, int maxiter, double tol)
{
    //std::cout << "Solver: maxiter = " << maxiter << ", tolerance = "
    //    << tol << std::endl;
    const int ncell = ebop.index.size();
    int iter = 0;
    double resid = tol;
    double param = .125; // =1/8, only valid for this problem

    //std::cout << "   beginning solver point Jacobi iterations ..." << std::endl;
    double maxresid;
    while (iter < maxiter)
    {
        zero(tmp);
        applybcsr(ebbcsr, soln, tmp);
        // Point Jacobi update
        maxresid = 0;
        for (int n = 0; n < ncell; ++n)
        {
            const size_t ix = ebop.index[n];
            const double resid = tmp.ptr[ix] - rhs.ptr[ix];
            maxresid = std::max<double>(maxresid, abs(resid));
            soln.ptr[ix] += resid * param;
        }

        if (maxresid < tol)
            break;
        ++iter;
    }
    //std::cout << "Solver exited, " << iter << " iterations, residual = "
    //    << maxresid << std::endl;
}

// Check if the sparse operator output matches the expected value 0
bool check(arr_t& exact, arr_t& soln, arr_t& opval, ebix_t& eb, 
    double tol, bool verbose)
{
  double truth = 4.*dx*dx; // Laplacian of (x^2 + y^2) * dx^2
  bool pass = true;
  for (int j=lo[1]; j <= hi[1]; ++j)
  for (int i=lo[0]; i <= hi[0]; ++i)
  {
    if (eb.inout[j][i] < 0)
      continue;
    double diff = opval.data[j][i] - truth;
    if (abs(diff) > tol)
    {
      if (verbose)
        std::cout << "Opval error in cell "
          << "[" << j << "][" << i << "]" << std::endl
          << "  --> expected " << truth 
          << ", value = " << opval.data[j][i] 
          << ", diff = " << abs(diff) << std::endl;
      pass = false;
    }
    diff = exact.data[j][i] - soln.data[j][i];
    if (abs(diff) > tol)
    {
      if (verbose)
        std::cout << "Solution error in cell "
          << "[" << j << "][" << i << "]" << std::endl
          << "  --> expected " << exact.data[j][i]
          << ", value = " << soln.data[j][i] 
          << ", diff = " << abs(diff) << std::endl;
      pass = false;
    }
  }
  return pass;
}

bool checkdiff(arr_t& val1, arr_t& val2, const double& tol)
{
  bool pass=true;
  for (size_t i=0; i < NTOT; ++i)
  {
    double diff = val1.ptr[i] - val2.ptr[i];
    if (abs(diff) > tol)
    {
      std::cout << "Diff of arrays in entry " << "[" << i 
        << "] greater than " << tol << ", diff = " << abs(diff) << std::endl;
      pass=false;
    }
  }
  return pass;
}

int main (int argc, char* argv[])
{
  arr_t exact, soln, opval, tmp, bcsrtmp;
  ebix_t eb;
  ebsparseop_t ebop;
  ebsparsematrix_t ebmat;
  ebbcsrmatrix_t ebbcsr;
  ebcsrmatrix_t ebcsr;
  double tol=1e-14; // floating point tolerance for tests

//test input
std::vector<size_t> rows{ 0,    0,   0,   0,   1,  1,    1,   1,   2,   2,   2,      4,  4, 5, 0, 63};//{ 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
std::vector<size_t> cols{ 0,    1,   2,   3,   0,  1,    2,   4,   0,   3,   4,      2,  5, 5, 63, 63};//{ 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
std::vector<double> vals{ 1,    2,   5,   6,   3,  4,    7,   8,   9,   1,   2,      4,  5, 6, 63, 2};//{0,0,2,5,0,0,21,22,1,1,1,1,1,1,1,1};

  ebsparsematrix_t mat;
  mat.N=64;
  mat.nzrow = rows;
  mat.nzcol = cols;
  mat.entry = vals;



  //initialize(exact, eb, ebop, ebmat);
  // dumpmat(ebmat, false);
  // dump(soln, std::string("solution value"));
  // dumpop(eb, ebop);
  convertToBCSR(mat, ebbcsr);//**ebmat

  //Compact values
  std::vector<double> multiplier(mat.N, 0);//**ebmat
  for(int i=0; i<mat.N; i++)//**ebmat
    {
	multiplier[i]=i;
    }
  ebbcsrBMmatrix_t bmp;
  generateBitmap(ebbcsr, bmp);
  std::vector<double> hOut;
  hOut = calculateTime(bmp, multiplier);
  /*std::vector<double> dOut;
  dOut = compactvaluescu(bmp, multiplier);
  if(hOut.size()!=dOut.size()){
    std::cout<<std::endl<<"Compressed Block: Incorect data. Host and Device data do not match.TEST FAILED!";
  }
  else{
  int ctr=0;
  for(int i=0; i<dOut.size(); i++)
     (dOut[i]==hOut[i]?ctr++:ctr);
  (ctr==dOut.size()?std::cout<<std::endl<<"Compressed Block: GPU and CPU data are consistent. TEST PASSED!"<<std::endl:std::cout<<std::endl<<"Compressed Block: GPU and CPU data are not consistent. TEST FAILED!");
  }*/
  for(int j=0; j<hOut.size(); j++){
    printf("%d %f %f\n",j, hOut[j],hOut[j]);}

  return 0; 
  // Check that the matrix and operator produce the same result
  zero(opval);
  applyop(ebop, exact, opval, false);
  zero(tmp);
  applymat(ebmat, exact, tmp);

  zero(bcsrtmp);
  applybcsr(ebbcsr, exact, bcsrtmp);


  //if (checkdiff(opval, tmp, tol))
  //  std::cout << "Matrix vs. operator test passed!" << std::endl;
  //else
  //  std::cout << "Matrix vs. operator test failed!" << std::endl;

  if (checkdiff(bcsrtmp, tmp, tol))
    std::cout << std::endl<<"Matrix vs. BCSR test passed!" << std::endl;
  else
    std::cout << "Matrix vs. BCSR test failed!" << std::endl;

  // Calculate the rhs for the solver, using the operator
  zero(opval);
  applyop(ebop, exact, opval, true);
  // Alternatively, use the eb operator
  // applyop(ebop, exact, opval, true);
  // dump(opval, std::string("operator value"));
  // zero(soln); // initial guess
  // copy(soln, exact); // use this to check for 0 iterations
  int maxiter = 20*NX*NX; // theory says it should scale like this
  // solve(ebop, soln, opval, tmp, maxiter, tol); // solve using the operator
  // dump(opval, std::string("opval after solve"));

  // Now try with the sparse matrix op - rhs has already been bc-fixed
  
  int num_runs = 100;
  int num_ignore = 10;
  TimeLog::ignore(num_ignore);
  for (int i = 0; i < num_runs; i++) {
      zero(soln); // initial guess
      zero(opval); // to zero out invalid region values 
      applymat(ebmat, exact, opval);
      auto t1 = std::chrono::high_resolution_clock::now();
      solvesparse(ebmat, ebop, soln, opval, tmp, maxiter, tol); // solve using spmv 
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> ms_sparse = t2 - t1;

      TimeLog::log(ms_sparse.count());
  }
  double sparseavg = TimeLog::calc_avg();
  TimeLog::reset();
 
  // dump(soln, std::string("solution value"));

  // Now try with the bcsr op - rhs has already been bc-fixed
  TimeLog::ignore(num_ignore);
  for (int i = 0; i < num_runs; i++) {
      zero(soln); // initial guess
      zero(opval); // to zero out invalid region values 
      applybcsr(ebbcsr, exact, opval);
      auto t1 = std::chrono::high_resolution_clock::now();
      solvebcsr(ebbcsr, ebop, soln, opval, tmp, maxiter, tol); // solve using spmv 
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> ms_bcsr = t2 - t1;
      TimeLog::log(ms_bcsr.count());
  }
  double bcsravg = TimeLog::calc_avg();
  // dump(soln, std::string("solution value"));

  // Check vs. expected result
  applyop(ebop, soln, opval, true);
  // applybc(ebop, opval, +1);
  tol*=NX*NX; // theory says it should scale like this
  if (check(exact, soln, opval, eb, tol, true))
    std::cout << "Test passed!" << std::endl;
  else
    std::cout << "Test failed!" << std::endl;

  std::cout << "Solve sparse time: " << sparseavg << "ms\n";
  std::cout << "Solve bcsr time: " << bcsravg << "ms\n";
  double percentdiff = ((sparseavg - bcsravg) / ((sparseavg + bcsravg) / 2)) * 100;
  std::cout << "Percentage difference of BCSR implementation: " << percentdiff << "%";

  return 0;
}
