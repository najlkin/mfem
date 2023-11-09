// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_BLOCKVECTOR
#define MFEM_BLOCKVECTOR

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "vector.hpp"

namespace mfem
{

//! @class BlockVector
/**
 * \brief A class to handle Vectors in a block fashion
 *
 * All data is contained in Vector::data, while blockVector is just a viewer for
 * this data.
 *
 */
class BlockVector: public Vector
{
protected:

   //! Number of blocks in the blockVector
   int numBlocks;
   //! Offset for each block start. (length numBlocks+1)
   /**
    * blockOffsets[i+1] - blockOffsets[i] is the size of block i.
    *
    * This array is not owned.
    */
   const int *blockOffsets;
   //! array of Vector objects used to extract blocks without allocating memory.
   /** This array is owned. */
   Vector *blocks;

   void SetBlocks();

public:
   //! empty constructor
   BlockVector();

   //! Constructor
   /**
    * bOffsets is an array of integers (length nBlocks+1) that tells the offsets
    * of each block start.
    */
   BlockVector(const Array<int> & bOffsets);

   /// Construct a BlockVector with the given MemoryType @a mt.
   BlockVector(const Array<int> & bOffsets, MemoryType mt);

   //! Copy constructor
   BlockVector(const BlockVector & block);

   //! View constructor
   /**
    * data is an array of double of length at least blockOffsets[numBlocks] that
    * contain all the values of the monolithic vector.  bOffsets is an array of
    * integers (length nBlocks+1) that tells the offsets of each block start.
    * nBlocks is the number of blocks.
    */
   BlockVector(double *data, const Array<int> & bOffsets);

   /// Wrap a Vector as a BlockVector
   BlockVector(Vector &v, const Array<int> &bOffsets);

   //! Return the number of blocks
   int NumBlocks() const { return numBlocks; }

   //! Assignment operator. this and original must have the same block structure.
   BlockVector & operator=(const BlockVector & original);
   //! Set each entry of this equal to val
   BlockVector & operator=(double val);

   //! Destructor
   ~BlockVector();

   //! Get the i-th vector in the block.
   Vector & GetBlock(int i) { return blocks[i]; }
   //! Get the i-th vector in the block (const version).
   const Vector & GetBlock(int i) const { return blocks[i]; }

   //! Get the i-th vector in the block
   void GetBlockView(int i, Vector & blockView);

   int BlockSize(int i) { return blockOffsets[i+1] - blockOffsets[i]; }

   //! Update method
   /**
    * data is an array of double of length at least blockOffsets[numBlocks] that
    * contain all the values of the monolithic vector.  bOffsets is an array of
    * integers (length nBlocks+1) that tells the offsets of each block start.
    * nBlocks is the number of blocks.
    */
   void Update(double *data, const Array<int> & bOffsets);

   void Update(Vector & data, const Array<int> & bOffsets);

   /// Update a BlockVector with new @a bOffsets and make sure it owns its data.
   /** The block-vector will be re-allocated if either:
       - the offsets @a bOffsets are different from the current offsets, or
       - currently, the block-vector does not own its data. */
   void Update(const Array<int> &bOffsets);

   /** @brief Update a BlockVector with new @a bOffsets and make sure it owns
       its data and uses the MemoryType @a mt. */
   /** The block-vector will be re-allocated if either:
       - the offsets @a bOffsets are different from the current offsets, or
       - currently, the block-vector does not own its data, or
       - currently, the block-vector does not use MemoryType @a mt. */
   void Update(const Array<int> &bOffsets, MemoryType mt);

   /** @brief Synchronize the memory location flags (i.e. the memory validity
       flags) of the big/monolithic block-vector with its sub-vector blocks. The
       big/monolithic vector has the correct memory location flags. */
   /** This method will copy the data validity flags from the big/monolithic
       block-vector to its sub-vector block. */
   void SyncToBlocks() const;

   /** @brief Synchronize the memory location flags (i.e. the memory validity
       flags) of the big/monolithic block-vector with its sub-vector blocks. The
       sub-vector blocks have the correct memory location flags. */
   /** This method will copy/move the data of the sub-vector blocks (if
       necessary) so that each block matches the memory location flags of the
       big/monolithic block-vector. */
   void SyncFromBlocks() const;
};

//! @class BlockRefVector
/**
 * \brief A class to handle Vectors in a block fashion using references
 */
class BlockRefVector: public AbstractVector
{
protected:

   //! Array of Vector objects
   Array<Vector*> blocks;

   //! Array of ownership flags
   Array<bool> block_own;

public:
   //! empty constructor
   BlockRefVector() { }

   //! Constructor
   BlockRefVector(const Array<Vector*> &blocks, bool bown=false);

   //! Constructor
   BlockRefVector(const Array<Vector*> &blocks, const Array<bool> &own);

   //! Constructor
   /**
    * bOffsets is an array of integers (length nBlocks+1) that tells the offsets
    * of each block start. Blocks are allocated and owned.
    */
   BlockRefVector(const Array<int> &bOffsets);

   //! Copy constructor (deep copy)
   BlockRefVector(const BlockRefVector & block);

   //! Return the number of blocks
   inline int NumBlocks() const { return blocks.Size(); }

   //! Set number of blocks
   /** @note The content is deleted when reallocating */
   void SetNumBlocks(int s) { if(s != blocks.Size()) { Destroy(); Init(s); } }

   //! Conversion to a monolithic vector
   operator Vector(); 

   //! Assignment operator
   /** This and original must have the same block structure. */
   BlockRefVector & operator=(const BlockRefVector & original);

   ///! Assignment opertor
   /** The monolithic vector is split and copied by parts. */
   BlockRefVector & operator=(const Vector &v);

   ///! Assignment opertor
   /** The monolithic data are split and copied by parts. */
   BlockRefVector & operator=(const double *v);

   BlockRefVector & operator=(double val)
   { for(Vector *v : blocks) if(v) *v = val; return *this; }

   BlockRefVector & operator+=(const BlockRefVector &x)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) *blocks[i] += x.GetBlock(i); return *this; }

   BlockRefVector & operator+=(double val)
   { for(Vector *v : blocks) if(v) *v += val; return *this; }

   BlockRefVector & operator-=(const BlockRefVector &x)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) *blocks[i] -= x.GetBlock(i); return *this; }
   
   BlockRefVector & operator-=(double val)
   { for(Vector *v : blocks) if(v) *v -= val; return *this; }
   
   BlockRefVector & operator*=(const BlockRefVector &x)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) *blocks[i] *= x.GetBlock(i); return *this; }

   BlockRefVector & operator*=(double val)
   { for(Vector *v : blocks) if(v) *v *= val; return *this; }
   
   BlockRefVector & operator/=(const BlockRefVector &x)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) *blocks[i] /= x.GetBlock(i); return *this; }
   
   BlockRefVector & operator/=(double val)
   { for(Vector *v : blocks) if(v) *v /= val; return *this; }

    /// (*this) += a * Va
   BlockRefVector &Add(const double a, const BlockRefVector &Va)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) blocks[i]->Add(a, Va.GetBlock(i)); return *this; }
   
   /// (*this) = a * x
   BlockRefVector &Set(const double a, const BlockRefVector &x)
   { for(int i = 0; i < blocks.Size(); i++) if(blocks[i]) blocks[i]->Set(a, x.GetBlock(i)); return *this; }

   /// (*this) = -(*this)
   void Neg() { for(Vector *v : blocks) if(v) v->Neg(); }

   //! Get the i-th vector in the block.
   Vector & GetBlock(int i) { MFEM_ASSERT(blocks[i],"") return *blocks[i]; }

   //! Get the i-th vector in the block (const version).
   const Vector & GetBlock(int i) const { MFEM_ASSERT(blocks[i],"") return *blocks[i]; }

   //! Set the i-th vector in the block.
   void SetBlock(int i, Vector *vec, bool bown=false);

   //! Get the i-th vector in the block
   void GetBlockView(int i, Vector & blockView)
   { if(blocks[i]) blockView.MakeRef(*blocks[i], 0, blocks[i]->Size()); }

   //! Check if the block is zero
   bool IsZeroBlock(int i) { return (blocks[i] == NULL); }

   //! Get size of the the i-th block
   int BlockSize(int i) { return (blocks[i])?(blocks[i]->Size()):(0); }

   void SetSize(int s) { MFEM_ABORT("Not implemented for this class");  }

   double& Elem(int i);
   const double& Elem(int i) const;

   //! Update method
   /**
    * data is an array of double of length at least blockOffsets[numBlocks] that
    * contain all the values of the monolithic vector.  bOffsets is an array of
    * integers (length nBlocks+1) that tells the offsets of each block start.
    * nBlocks is the number of blocks. Owned vectors are made, pointing to
    * the data.
    */
   void Update(double *data, const Array<int> & bOffsets);
   void Update(Vector &data, const Array<int> & bOffsets);

   //! Update method (owned vectors)
   /** The same as Update(double*, const Array<int> &), but only owned vectors
    * are replaced by those pointing to the provided data. */
   void UpdateOwned(double *data, const Array<int> & bOffsets);
   void UpdateOwned(Vector &data, const Array<int> & bOffsets);

   //! Update method (reference)
   /** Makes a complete reference (no ownership) to the given vector. */
   void Update(BlockRefVector &data);

   //! Update a BlockVector with new @a bOffsets and make sure it owns its data.
   /** The block-vector will be re-allocated if either:
       - the offsets @a bOffsets are different from the current offsets, or
       - currently, the block-vector does not own its data. */
   void Update(const Array<int> &bOffsets);

   //! Update a BlockVector with new @a bOffsets for the owned data only.
   /** @see Update(const Array<int> &) */
   void UpdateOwned(const Array<int> &bOffsets);

   //! Destroy the vectors
   void Destroy();

   //! Destructor
   ~BlockRefVector() { Destroy(); }

private:
   void Init(int s);
   void RecalcSize();
};

inline void add(const BlockRefVector &v1, const BlockRefVector &v2, BlockRefVector &v)
{
   for(int i = 0; i < v.NumBlocks(); i++)
      if(!v.IsZeroBlock(i))
         add(v1.GetBlock(i), v2.GetBlock(i), v.GetBlock(i));
}

inline void add(const BlockRefVector &v1, double alpha, const BlockRefVector &v2, BlockRefVector &v)
{
   for(int i = 0; i < v.NumBlocks(); i++)
      if(!v.IsZeroBlock(i))
         add(v1.GetBlock(i), alpha, v2.GetBlock(i), v.GetBlock(i));
}

inline void add(const double a, const BlockRefVector &x, const BlockRefVector &y, BlockRefVector &z)
{
   for(int i = 0; i < z.NumBlocks(); i++)
      if(!z.IsZeroBlock(i))
         add(a, x.GetBlock(i), y.GetBlock(i), z.GetBlock(i));
}

inline void add(const double a, const BlockRefVector &x, const double b, const BlockRefVector &y, BlockRefVector &z)
{
   for(int i = 0; i < z.NumBlocks(); i++)
      if(!z.IsZeroBlock(i))
         add(a, x.GetBlock(i), b, y.GetBlock(i), z.GetBlock(i));
}

inline void subtract(const BlockRefVector &v1, const BlockRefVector &v2, BlockRefVector &v)
{
   for(int i = 0; i < v.NumBlocks(); i++)
      if(!v.IsZeroBlock(i))
         subtract(v1.GetBlock(i), v2.GetBlock(i), v.GetBlock(i));
}

inline void subtract(const double a, const BlockRefVector &x, const BlockRefVector &y, BlockRefVector &z)
{
   for(int i = 0; i < z.NumBlocks(); i++)
      if(!z.IsZeroBlock(i))
         subtract(a, x.GetBlock(i), y.GetBlock(i), z.GetBlock(i));
}

}

#endif /* MFEM_BLOCKVECTOR */
