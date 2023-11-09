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

#include "../general/array.hpp"
#include "vector.hpp"
#include "blockvector.hpp"

namespace mfem
{

void BlockVector::SetBlocks()
{
   for (int i = 0; i < numBlocks; ++i)
   {
      blocks[i].MakeRef(*this, blockOffsets[i], BlockSize(i));
   }
}

BlockVector::BlockVector():
   Vector(),
   numBlocks(0),
   blockOffsets(NULL),
   blocks(NULL)
{

}

//! Standard constructor
BlockVector::BlockVector(const Array<int> & bOffsets):
   Vector(bOffsets.Last()),
   numBlocks(bOffsets.Size()-1),
   blockOffsets(bOffsets.GetData())
{
   blocks = new Vector[numBlocks];
   SetBlocks();
}

BlockVector::BlockVector(const Array<int> & bOffsets, MemoryType mt)
   : Vector(bOffsets.Last(), mt),
     numBlocks(bOffsets.Size()-1),
     blockOffsets(bOffsets.GetData())
{
   blocks = new Vector[numBlocks];
   SetBlocks();
}

//! Copy constructor
BlockVector::BlockVector(const BlockVector & v):
   Vector(v),
   numBlocks(v.numBlocks),
   blockOffsets(v.blockOffsets)
{
   blocks = new Vector[numBlocks];
   SetBlocks();
}

//! View constructor
BlockVector::BlockVector(double *data, const Array<int> & bOffsets):
   Vector(data, bOffsets.Last()),
   numBlocks(bOffsets.Size()-1),
   blockOffsets(bOffsets.GetData())
{
   blocks = new Vector[numBlocks];
   SetBlocks();
}

BlockVector::BlockVector(Vector &v, const Array<int> &bOffsets)
   : Vector(),
     numBlocks(bOffsets.Size()-1),
     blockOffsets(bOffsets.GetData())
{
   MakeRef(v, 0, blockOffsets[numBlocks]);
   blocks = new Vector[numBlocks];
   SetBlocks();
}

void BlockVector::Update(double *data, const Array<int> & bOffsets)
{
   NewDataAndSize(data, bOffsets.Last());
   blockOffsets = bOffsets.GetData();
   if (numBlocks != bOffsets.Size()-1)
   {
      delete [] blocks;
      numBlocks = bOffsets.Size()-1;
      blocks = new Vector[numBlocks];
   }
   SetBlocks();
}

void BlockVector::Update(Vector & data, const Array<int> & bOffsets)
{
   blockOffsets = bOffsets.GetData();
   if (numBlocks != bOffsets.Size()-1)
   {
      delete [] blocks;
      numBlocks = bOffsets.Size()-1;
      blocks = new Vector[numBlocks];
   }

   for (int i = 0; i < numBlocks; ++i)
   {
      blocks[i].MakeRef(data, blockOffsets[i], BlockSize(i));
   }
   MakeRef(data, 0, blockOffsets[numBlocks]);
}

void BlockVector::Update(const Array<int> &bOffsets)
{
   Update(bOffsets, data.GetMemoryType());
}

void BlockVector::Update(const Array<int> &bOffsets, MemoryType mt)
{
   blockOffsets = bOffsets.GetData();
   if (OwnsData() && data.GetMemoryType() == mt)
   {
      // check if 'bOffsets' agree with the 'blocks'
      if (bOffsets.Size() == numBlocks+1)
      {
         if (numBlocks == 0) { return; }
         if (Size() == bOffsets.Last())
         {
            for (int i = numBlocks - 1; true; i--)
            {
               if (i < 0) { return; }
               if (blocks[i].Size() != bOffsets[i+1] - bOffsets[i]) { break; }
               MFEM_ASSERT(blocks[i].GetData() == data + bOffsets[i],
                           "invalid blocks[" << i << ']');
            }
         }
      }
   }
   else
   {
      Destroy();
   }
   SetSize(bOffsets.Last(), mt);
   if (numBlocks != bOffsets.Size()-1)
   {
      delete [] blocks;
      numBlocks = bOffsets.Size()-1;
      blocks = new Vector[numBlocks];
   }
   SetBlocks();
}

BlockVector & BlockVector::operator=(const BlockVector & original)
{
   if (numBlocks!=original.numBlocks)
   {
      mfem_error("Number of Blocks don't match in BlockVector::operator=");
   }

   for (int i(0); i <= numBlocks; ++i)
   {
      if (blockOffsets[i]!=original.blockOffsets[i])
      {
         mfem_error("Size of Blocks don't match in BlockVector::operator=");
      }
   }

   Vector::operator=(original);

   return *this;
}

BlockVector & BlockVector::operator=(double val)
{
   Vector::operator=(val);
   return *this;
}

//! Destructor
BlockVector::~BlockVector()
{
   delete [] blocks;
}

void BlockVector::GetBlockView(int i, Vector & blockView)
{
   blockView.MakeRef(*this, blockOffsets[i], BlockSize(i));
}

void BlockVector::SyncToBlocks() const
{
   for (int i = 0; i < numBlocks; ++i)
   {
      blocks[i].SyncMemory(*this);
   }
}

void BlockVector::SyncFromBlocks() const
{
   for (int i = 0; i < numBlocks; ++i)
   {
      blocks[i].SyncAliasMemory(*this);
   }
}

BlockRefVector::BlockRefVector(const Array<Vector*> &blocks_, bool bown)
{
   blocks = blocks_;
   block_own = bown;
   RecalcSize();
}

BlockRefVector::BlockRefVector(const Array<Vector*> &blocks_, const Array<bool> &owns)
{
   blocks = blocks_;
   block_own = owns;
   RecalcSize();
}

BlockRefVector::BlockRefVector(const Array<int> &offsets)
{
   Update(offsets);
}

BlockRefVector::BlockRefVector(const BlockRefVector &orig)
{
   Init(orig.NumBlocks());
   for(int i = 0; i < blocks.Size(); i++)
      if(orig.blocks[i])
      {
         blocks[i] = new Vector(*orig.blocks[i]);
         block_own[i] = true;
      }
   size = orig.size;
}

BlockRefVector::operator Vector()
{
   Vector v(Size());
   Vector view;
   int ioff = 0;
   for(Vector *b : blocks)
   {
      if(!b)
         continue;
      view.MakeRef(v, ioff, b->Size());
      view = *b;
      ioff += b->Size();
   }
   return v;
}

BlockRefVector &BlockRefVector::operator=(const BlockRefVector &orig)
{
   for(int i = 0; i < blocks.Size(); i++)
      if(orig.blocks[i])
         *blocks[i] = *orig.blocks[i];

   return *this;
}

BlockRefVector &BlockRefVector::operator=(const Vector &v)
{
   MFEM_ASSERT(v.Size() >= Size(), "Insufficiently long vector");
   return operator=(v.GetData());
}

BlockRefVector &BlockRefVector::operator=(const double *v)
{
   int ioff = 0;
   for(Vector *b : blocks)
   {
      if(!b)
         continue;
      *b = v + ioff;
      ioff += b->Size();
   }

   return *this;
}

void BlockRefVector::SetBlock(int i, Vector *vec, bool bown)
{
   if(blocks[i] == vec)
   {
      block_own[i] = bown;
      return;
   }

   if(blocks[i])
   {
      size -= blocks[i]->Size();
      if(block_own[i])
         delete blocks[i];
   }
   blocks[i] = vec;
   block_own[i] = bown;
   if(vec)
      size += vec->Size();
}

double& BlockRefVector::Elem(int i)
{
   int b = 0;
   int ioff = 0;
   for(; b < blocks.Size(); b++)
   {
      if(!blocks[b])
         continue;
      if(i < ioff + blocks[b]->Size())
         break;
      ioff += blocks[b]->Size();
   }
   return blocks[b]->Elem(i-ioff);
}

const double& BlockRefVector::Elem(int i) const
{
   int b = 0;
   int ioff = 0;
   for(; b < blocks.Size(); b++)
   {
      if(!blocks[b])
         continue;
      if(i < ioff + blocks[b]->Size())
         break;
      ioff += blocks[b]->Size();
   }
   return blocks[b]->Elem(i-ioff);
}

void BlockRefVector::Update(double *data, const Array<int> &bOffsets)
{
   if(blocks.Size() != bOffsets.Size()-1)
   {
      Destroy();
      Init(bOffsets.Size()-1);
   }

   for(int i = 0; i < bOffsets.Size()-1; i++)
   {
      const int size = bOffsets[i+1] - bOffsets[i];
      if(blocks[i] && block_own[i])
      {
         if(blocks[i]->Size() == size
         && blocks[i]->GetData() == data + bOffsets[i])
            continue;
      }
      if(block_own[i])
         delete blocks[i];
      if(size > 0)
         blocks[i] = new Vector(data + bOffsets[i], size);
      else
         blocks[i] = NULL;
      block_own[i] = true;
   }

   size = bOffsets.Last();
}


void BlockRefVector::Update(Vector &data, const Array<int> &bOffsets)
{
   MFEM_ASSERT(data.Size() == bOffsets.Last(), "");
   Update(data.GetData(), bOffsets);
}

void BlockRefVector::UpdateOwned(double *data, const Array<int> &bOffsets)
{
   int nowned = 0;
   for(bool own : block_own) { if(own) nowned++; }

   if(nowned != bOffsets.Size()-1)
   {
      Destroy();
      Init(bOffsets.Size()-1);
      block_own = true;
   }

   int ioff = 0;
   bool recalc_size = false;
   for(int i = 0; i < blocks.Size(); i++)
   {
      if(!block_own[i])
         continue;
      const int size = bOffsets[ioff+1] - bOffsets[ioff];
      ioff++;
      if(blocks[i])
      {
         if(blocks[i]->Size() == size
         && blocks[i]->GetData() == data + bOffsets[ioff-1])
            continue;
      }
      delete blocks[i];
      if(size > 0)
         blocks[i] = new Vector(data + bOffsets[ioff-1], size);
      else
         blocks[i] = NULL;
      recalc_size = true;
   }

   if(recalc_size)
      RecalcSize();
}

void BlockRefVector::UpdateOwned(Vector &data, const Array<int> &bOffsets)
{
   MFEM_ASSERT(data.Size() == bOffsets.Last(), "");
   UpdateOwned(data.GetData(), bOffsets);
}

void BlockRefVector::Update(BlockRefVector &data)
{
   if(data.NumBlocks() != blocks.Size())
   {
      Destroy();
      Init(data.NumBlocks());
   }

   for(int i = 0; i < data.NumBlocks(); i++)
   {
      const int size = data.BlockSize(i);
      if(blocks[i])
      {
         if(blocks[i]->Size() == size
         && blocks[i] == &data.GetBlock(i))
         {
            block_own[i] = false;
            continue;
         }
      }
      if(block_own[i])
         delete blocks[i];
      if(size > 0)
         blocks[i] = &data.GetBlock(i);
      else
         blocks[i] = NULL;
      block_own[i] = false;
   }

   size = data.Size();
}

void BlockRefVector::Update(const Array<int> &bOffsets)
{
   if(blocks.Size() != bOffsets.Size()-1)
   {
      Destroy();
      Init(bOffsets.Size()-1);
   }

   for(int i = 0; i < bOffsets.Size()-1; i++)
   {
      const int size = bOffsets[i+1] - bOffsets[i];
      if(blocks[i] && block_own[i])
      {
         if(blocks[i]->Size() == size)
            continue;
      }
      if(block_own[i])
         delete blocks[i];
      if(size > 0)
         blocks[i] = new Vector(size);
      else
         blocks[i] = NULL;
      block_own[i] = true;
   }

   size = bOffsets.Last();
}

void BlockRefVector::UpdateOwned(const Array<int> &bOffsets)
{
   int nowned = 0;
   for(bool own : block_own) { if(own) nowned++; }
   
   if(nowned != bOffsets.Size()-1)
   {
      Destroy();
      Init(bOffsets.Size()-1);
      block_own = true;
   }

   int ioff = 0;
   bool recalc_size = false;
   for(int i = 0; i < blocks.Size(); i++)
   {
      if(!block_own[i])
         continue;
      const int size = bOffsets[ioff+1] - bOffsets[ioff];
      ioff++;
      if(blocks[i])
      {
         if(blocks[i]->Size() == size)
            continue;
      }
      delete blocks[i];
      if(size > 0)
         blocks[i] = new Vector(size);
      else
         blocks[i] = NULL;
      recalc_size = true;
   }

   if(recalc_size)
      RecalcSize();
}

void BlockRefVector::Destroy()
{
   for(int i = 0; i < blocks.Size(); i++)
   {
      if(block_own[i])
         delete blocks[i];
   }
   blocks.SetSize(0);
   block_own.SetSize(0);
   size = 0;
}

void BlockRefVector::Init(int s)
{
   blocks.SetSize(s);
   blocks = NULL;
   block_own.SetSize(s);
   block_own = false;
   size = 0;
}

void BlockRefVector::RecalcSize()
{
   size = 0;
   for(Vector *v : blocks)
      if(v)
         size += v->Size();
}

}
