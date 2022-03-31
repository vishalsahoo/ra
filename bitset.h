#ifndef BITSET_H
#define BITSET_H
#include "forcuda.h"
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <tuple>

template <typename T> struct BitSetRow{
  int firstbit, lastbit, datidx, rowinblock;
  T* bitset;
  typedef int IterValType;
  typedef BitSetRow<T> MyType;
  class iterator
      : public std::iterator<std::input_iterator_tag, IterValType, long,
                             const IterValType *, IterValType> {
    int cnt = -1;
    uint32_t val, val1, poscopy; bool retflag = false, val_change = true;
    int rowinblock_copy, bs_copy, elempos_start, elempos_end, ctr=-1;
    int pos=-1, ctz, elempos = -1, i=0, pos1=-1;

    BitSetRow<T> *container;

    void next() {
      if (container->bitset->bs < 32) {
          if (val==0 || i>=container->bitset->bs || retflag == true){
          pos = -1; poscopy = -1;
          return;
          }
          else{
            ctz = __builtin_ffs(val)-1;
            if (ctz == 0){
              cnt++; pos++; i++; val = val>>1; poscopy++;
              }
            else if (ctz == -1){
              retflag = true;
              } 
            else{
              val = val>>ctz;
              i += ctz; pos += ctz; poscopy += ctz; next();
            }
          }
      }
      else {
        if(val_change){
          elempos_start += 1; ctr += 1;
          val1 = container->bitset->dat[elempos_start]; retflag=false; val_change=false; pos1 = -1;
          }
        if (elempos_start>=elempos_end){
          pos1 = -1; poscopy = -1; pos=pos1;
          return;
          }
        else{
            ctz = __builtin_ffs(val1)-1;
            if (ctz == 0){
              cnt++; pos1++; i++; val1 = val1>>1; poscopy++; pos++;
              }
            else if (ctz == -1){
              val_change = true; next();
              } 
            else{
              val1 = val1>>ctz;
              i += ctz; pos1 += ctz; poscopy += ctz; 
              pos=(32*ctr)+pos1; 
              next(); 
            }
          }
        }
    }

  public:
    explicit iterator(MyType *container, bool end = false)
        : container(container) {
      container->firstbit = (container->bitset->bs * container->rowinblock)%container->bitset->vlen;
      container->lastbit = container->firstbit + (container->bitset->bs-1);
      uint32_t &s = container->bitset->dat[(container->bitset->bs*container->rowinblock)/container->bitset->vlen];
      uint32_t mask = ~(~0 << (container->lastbit - container->firstbit + 1));
      val = (s >> container->firstbit) & mask;
      rowinblock_copy = container->rowinblock;
      elempos_start = (container->rowinblock*((container->bitset->bs)/32))-1;
      elempos_end = container->rowinblock*((container->bitset->bs)/32) + ((container->bitset->bs)/32);
      next();
    }

    iterator& operator++() {
      next(); return *this;
    }

    bool operator==(iterator other) const {
      return 0;
    }

    bool operator!=(iterator other) const {
     return (poscopy!=-1);
    }

    std::tuple<int, int> operator*() const {
      return std::tuple<int, int> (cnt, pos);
    }
  };

  iterator begin() {return iterator(this); }
  iterator end() {return iterator(this, false); }
};

template <unsigned rows, unsigned cols> struct BitSet {
  static constexpr unsigned vlen = 32;
  static constexpr unsigned sz = rows * cols;
  static constexpr unsigned len = (sz + vlen - 1) / vlen;
  static constexpr unsigned r = rows;
  static constexpr unsigned c = cols;
  static constexpr unsigned bs = cols;
  typedef BitSet<rows, cols> MyType;
  uint32_t dat[len];

  BitSetRow<MyType> operator [](int row) {
    BitSetRow<MyType> obj1;
    obj1.datidx = (bs*row)/32;
    obj1.firstbit = bs*row;
    obj1.lastbit = bs*(row+1)-1;
    obj1.bitset = this;
    obj1.rowinblock = row;
    return obj1;
  }

  void set(unsigned pos) {
    auto &s = dat[pos / vlen];
    auto p = 1 << (pos % vlen);
    s = s | p;
  }
#if defined(__CUDACC__)
  __device__ void set_atomic(unsigned pos) {
    auto p = 1 << (pos % vlen);
    atomicOr(&dat[pos / vlen], p);
  }
#endif
  void unset(unsigned pos) {
    auto &s = dat[pos / vlen];
    uint32_t p = 0;
    p = (~p) - (1 << (pos % vlen));
    s = s & p;
  }
#if defined(__CUDACC__)
  __device__ void unset_atomic(unsigned pos) {
    uint32_t p = 0;
    p = (~p) - (1 << (pos % vlen));
    atomicAnd(&dat[pos / vlen], p);
  }
#endif
  FORCUDA
  bool test(unsigned pos) {
    auto &s = dat[pos / vlen];
    auto p = 1 << (pos % vlen);
    return (s & p) > 0;
  }

#if defined(__CUDACC__)
  __device__ bool test_atomic(unsigned pos) {
    auto &s = dat[pos / vlen];
    auto p = 1 << (pos % vlen);
    return (s & p) > 0;
  }
#endif

  // CPU only
  BitSet() {
    for (int i = 0; i < (sz + vlen - 1) / vlen; ++i)
      dat[i] = 0;
  }

  std::vector<int> ctzCalc(int bRow, int bs) {
    // bs = BLOCKSIZE;
    std::vector<int> vectPos;
    if (bs <= vlen) {
      int datIdx = (bs * bRow) / vlen, datStart = (bs * bRow) % vlen,
          datEnd = datStart + (bs - 1);
      vectPos = opsBelow32(datIdx, datStart, datEnd, vectPos, 0, bs);
    } else {
      for (int ii = (bs / vlen) * bRow, whC = 0;
           ii < ((bs / vlen) * bRow) + (bs / vlen); ii++, whC++)
        vectPos = ops(ii, 0, vlen, vectPos, ii, bs, whC);
    }
    //	for(int xxx=0; xxx<vectPos.size(); xxx++) std::cout<<vectPos[xxx]<<" ";
    //	std::cout<<"\n";
    return vectPos;
  }

  std::vector<int> opsBelow32(int datIdx, int start, int end,
                              std::vector<int> vectPos, int off, int bs) {
    // int bs = BLOCKSIZE;
    auto &s = dat[datIdx];
    auto mask = ~(~0 << (end - start + 1));
    auto value = (s >> start) & mask;
    int i = 0;
    while (i < bs) {
      int ctz = __builtin_ctz(value);
      if (ctz == 0) {
        vectPos.push_back(off + i);
        // vectPos.push_back((whC*vlen)+i);
        value = value >> 1;
        i++;
      } else {
        value = value >> ctz;
        i += ctz;
      }
    }
    return vectPos;
  }

  std::vector<int> ops(int datIdx, int start, int end, std::vector<int> vectPos,
                       int off, int bs, int whC) {
    auto value = dat[datIdx];
    int i = 0;
    while (i < vlen) {
      int ctz = __builtin_ctz(value);
      if (ctz == 0) {
        vectPos.push_back((whC * vlen) + i);
        value = value >> 1;
        i++;
      } else {
        value = value >> ctz;
        i += ctz;
      }
    }
    return vectPos;
  }
};

#endif // BITSET_H