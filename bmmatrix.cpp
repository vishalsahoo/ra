//#include "common.h"
#include <iostream>
//#include "bitset.h"
//#include "bittemplate.h"
#include <bits/stdc++.h>
#include <tuple>
//#include <openacc.h>
using namespace std;

std::vector<double> calculateTime(ebbcsrBMmatrix_t bmp, std::vector<double> multiplier)
{
  int _size = 100;
  std::vector<double> y((bmp.block_row_ptr.size() - 1)*BLOCKSIZE, 0);
  int bs = BLOCKSIZE;

  double *values_ptr = &bmp.values[0];
  size_t *value_pos_ptr = &bmp.value_pos[0];
  size_t *cols_ptr = &bmp.cols[0];
  double *multiplier_ptr = &multiplier[0];
  double *y_ptr = &y[0];
  size_t *block_row_ptr_ptr = &bmp.block_row_ptr[0];
  BitSet<64, 64> *bit_ptr = &bmp.bit[0];

  #pragma acc parallel loop copyin(values_ptr[0:bmp.values.size()], value_pos_ptr[0:bmp.value_pos.size()], cols_ptr[0:bmp.cols.size()], multiplier_ptr[0:multiplier.size()], block_row_ptr_ptr[0:bmp.block_row_ptr.size()], bit_ptr[0:bmp.bit.size()]) copyout(y_ptr[0:y.size()])
  for(int idx=0; idx<bs*(bmp.block_row_ptr.size()-1); idx++){
       int row = idx % bs;
       int block_row = idx / bs;
       int first_block = block_row_ptr_ptr[block_row];
       int last_block = block_row_ptr_ptr[block_row + 1];
       y_ptr[idx]=0;
       for(int block=first_block; block < last_block; block++) {
        for(auto it: bit_ptr[block][row]){
          y_ptr[idx] += values_ptr[value_pos_ptr[block*bs+row]+get<0>(it)] * multiplier_ptr[cols_ptr[block]+ get<1>(it)];
        }   
       }
    }
return y;
}

void generateBitmap(ebbcsrmatrix_t &ebbcsr, ebbcsrBMmatrix_t &bmp)
{
  int flag_cont = 0, bitmap_idx=0, flat_bitmap_idx=0, blocksize = ebbcsr.blocksize;
  BitSet<BLOCKSIZE, BLOCKSIZE> temp;
  for (int i=0; i<ebbcsr.values.size(); i=i+(blocksize*blocksize)) 
  {
    int flag = 0; //std::map<int, int> mapIdxVal_temp;
    for(int j=i; j<(i+(blocksize*blocksize)); j++)
    {
    if(ebbcsr.values[j]==0)
      {
	temp.unset(j%(BLOCKSIZE*BLOCKSIZE));
      }
    else
      {
	temp.set(j%(BLOCKSIZE*BLOCKSIZE));
	flag = 1;
	bmp.values.push_back(ebbcsr.values[j]);
      }
   } 
  bmp.bit.push_back(temp);
  }

  for (int i=0, pos_ctr=0; i<ebbcsr.values.size(); i=i+(blocksize))
  {
    int count=0, flag_pos = 0;
    for(int j=i; j<i+blocksize; j++)
    {
      if(ebbcsr.values[j]!=0)
	{
	  pos_ctr++;
      if(flag_pos==0)
	{
	  bmp.value_pos.push_back(pos_ctr-1); flag_pos = 1;
	}
	
        }
      if(j==(i+blocksize-1) and flag_pos == 0)
	{
	  bmp.value_pos.push_back(-1);
	}
      
    }
  } 

  for(int i=0; i<ebbcsr.cols.size(); i++)
  {
    bmp.cols.push_back(ebbcsr.cols[i]*BLOCKSIZE);
  }
  for(int i=0; i<ebbcsr.block_row_ptr.size(); i++)
  {
    bmp.block_row_ptr.push_back(ebbcsr.block_row_ptr[i]);
  }
}