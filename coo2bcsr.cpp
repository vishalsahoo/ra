#include "common.h"

// Converts a COO matrix to a BCSR matrix
void convertToBCSR(ebsparsematrix_t &ebmat, ebbcsrmatrix_t &ebbcsr) {
  std::map<std::pair<int, int>, std::vector<double>> blockmap;

  for (int n = 0; n < ebmat.entry.size(); ++n) {
    // Retrieve information from COO matrix
    const int i = ebmat.nzrow[n];
    const int j = ebmat.nzcol[n];
    const double e = ebmat.entry[n];

    // Calculate block starting point
    const int ib = i / ebbcsr.blocksize;
    const int jb = j / ebbcsr.blocksize;
    std::pair<int, int> blockstart = std::pair<int, int>(ib, jb);

    // Calculate where the nz should be inside the block
    const int ii = i % ebbcsr.blocksize;
    const int jj = j % ebbcsr.blocksize;

    if (blockmap.find(blockstart) != blockmap.end()) {
      blockmap.at(blockstart)[ii * ebbcsr.blocksize + jj] = e;
    } else {
      std::vector<double> newvector = std::vector<double>(ebbcsr.blocksize * ebbcsr.blocksize, 0);

      newvector[ii * ebbcsr.blocksize + jj] = e;
      blockmap.insert({blockstart, newvector});
    }
  }
  std::map<std::pair<int, int>, std::vector<double>>::iterator it;

  // Pre allocate the space in the vectors
  ebbcsr.cols.reserve(blockmap.size());
  ebbcsr.values.reserve(blockmap.size() * ebbcsr.blocksize * ebbcsr.blocksize);
  ebbcsr.block_row_ptr = std::vector<size_t>(ebmat.N/ebbcsr.blocksize + 1, 0);

  // Getting the first non-zero row
  int prev_block_id_row = blockmap.begin()->first.first;
  int count = 0;
  for (it = blockmap.begin(); it != blockmap.end(); it++) {
    ebbcsr.values.insert(ebbcsr.values.end(), it->second.begin(), it->second.end());
    ebbcsr.cols.push_back(it->first.second); // Builds the column vector
    // Sets up the compressed row by counting the number of blocks in that row
    // Adds it to the row vector once it moves onto the next row.
    if (it->first.first != prev_block_id_row) {
      ebbcsr.block_row_ptr[prev_block_id_row + 1] = count;
      prev_block_id_row = it->first.first;
      count = 0;
    }
    count++;
  }
  // This is needed to add the last non-zero row
  ebbcsr.block_row_ptr[prev_block_id_row + 1] = count;

  // This loop finishes the compressed row format by creating the ranges
  // of where a row starts and ends in terms of blocks. e.g. the number of
  // block in row x should be row[x + 1] - rox[x]
  for (int i = 0; i < ebmat.N / ebbcsr.blocksize; i++) {
      ebbcsr.block_row_ptr[i + 1] += ebbcsr.block_row_ptr[i];
  }

}

// Converts a COO to CSR format. Used to test correctness of the compressed 
// sparse row of a 1x1 BCSR as they should be equivalent
void cootocsr(ebsparsematrix_t& ebmat, ebcsrmatrix_t& ebcsr) {
    ebcsr.cols = ebmat.nzcol;
    ebcsr.values = ebmat.entry;
    // Need N to know the original number of rows.
    ebcsr.row_ptr = std::vector<size_t>(ebmat.N + 1, 0);
    for (int i = 0; i < ebmat.entry.size(); i++) {
        ebcsr.row_ptr[ebmat.nzrow[i] + 1]++;
    }
    for (int i = 0; i < ebmat.N; i++) {
        ebcsr.row_ptr[i + 1] += ebcsr.row_ptr[i];
    }

}
