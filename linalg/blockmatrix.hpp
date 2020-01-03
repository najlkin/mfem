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

#ifndef MFEM_BLOCKMATRIX
#define MFEM_BLOCKMATRIX

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "../general/globals.hpp"
#include "vector.hpp"
#include "sparsemat.hpp"

namespace mfem
{

class BlockMatrix : public AbstractSparseMatrix
{
public:
   //! Constructor for square block matrices
   /**
     @param offsets  offsets that mark the start of each row/column block (size
                     nRowBlocks+1).
     @note BlockMatrix will not own/copy the data contained in offsets.
    */
   BlockMatrix(const Array<int> & offsets);
   //! Constructor for rectangular block matrices
   /**
     @param row_offsets  offsets that mark the start of each row block (size
                         nRowBlocks+1).
     @param col_offsets  offsets that mark the start of each column block (size
                         nColBlocks+1).
     @note BlockMatrix will not own/copy the data contained in offsets.
    */
   BlockMatrix(const Array<int> & row_offsets, const Array<int> & col_offsets);
   //! Set A(i,j) = mat
   void SetBlock(int i, int j, SparseMatrix * mat);
   //! Return the number of row blocks
   int NumRowBlocks() const {return nRowBlocks; }
   //! Return the number of column blocks
   int NumColBlocks() const {return nColBlocks; }
   //! Return a reference to block (i,j). Reference may be invalid if Aij(i,j) == NULL
   SparseMatrix & GetBlock(int i, int j);
   //! Return a reference to block (i,j). Reference may be invalid if Aij(i,j)
   //! == NULL. (const version)
   const SparseMatrix & GetBlock(int i, int j) const;
   //! Check if block (i,j) is a zero block.
   int IsZeroBlock(int i, int j) const {return (Aij(i,j)==NULL) ? 1 : 0; }
   //! Return the row offsets for block starts
   Array<int> & RowOffsets() { return row_offsets; }
   //! Return the columns offsets for block starts
   Array<int> & ColOffsets() { return col_offsets; }
   //! Return the row offsets for block starts (const version)
   const Array<int> & RowOffsets() const { return row_offsets; }
   //! Return the row offsets for block starts (const version)
   const Array<int> & ColOffsets() const { return col_offsets; }
   //! Return the number of non zeros in row i
   int RowSize(const int i) const;

   /// Eliminate the row and column @a rc from the matrix.
   /** Eliminates the column and row @a rc, replacing the element (rc,rc) with
       1.0. Assumes that element (i,rc) is assembled if and only if the element
       (rc,i) is assembled. If @a dpolicy is specified, the element (rc,rc) is
       treated according to that policy. */
   void EliminateRowCol(int rc, DiagonalPolicy dpolicy = DIAG_ONE);

   /** @brief Similar to EliminateRowCol(int, DiagonalPolicy) + save the
       eliminated entries into @a Ae so that (*this) + Ae is equal to the
       original matrix */
   void EliminateRowCol(int rc, BlockMatrix &Ae,
                        DiagonalPolicy dpolicy = DIAG_ONE);

   //! Symmetric elimination of the marked degree of freedom.
   /**
     @param ess_bc_dofs  marker of the degree of freedom to be eliminated
                         dof i is eliminated if @a ess_bc_dofs[i] = 1.
     @param sol          vector that stores the values of the degree of freedom
                         that need to be eliminated
     @param rhs          vector that stores the rhs of the system.
   */
   void EliminateRowCol(Array<int> & ess_bc_dofs, Vector & sol, Vector & rhs);

   //! Symmetric elimination of the marked degree of freedom.
   /**
     @param ess_bc_dofs  marker of the degree of freedom to be eliminated
                         dof i is eliminated if @a ess_bc_dofs[i] = 1.
     @param Ae           block matrix storing eliminated dofs with the same
                         structure as (*this)
     @param dpolicy      diagonal entries policy.
   */
   void EliminateRowCol(Array<int> & ess_bc_dofs, BlockMatrix &Ae,
                        DiagonalPolicy dpolicy = DIAG_ONE);

   ///  Finalize all the submatrices
   virtual void Finalize(int skip_zeros = 1) { Finalize(skip_zeros, false); }
   /// A slightly more general version of the Finalize(int) method.
   void Finalize(int skip_zeros, bool fix_empty_rows);

   //! Returns a monolithic CSR matrix that represents this operator.
   SparseMatrix * CreateMonolithic() const;
   //! Export the monolithic matrix to file.
   void PrintMatlab(std::ostream & os = mfem::out) const;

   /// @name Matrix interface
   ///@{

   /// Returns reference to a_{ij}.
   virtual double& Elem (int i, int j);
   /// Returns constant reference to a_{ij}.
   virtual const double& Elem (int i, int j) const;
   /// Returns a pointer to (approximation) of the matrix inverse.
   virtual MatrixInverse * Inverse() const
   {
      mfem_error("BlockMatrix::Inverse not implemented \n");
      return static_cast<MatrixInverse*>(NULL);
   }
   ///@}

   ///@name AbstractSparseMatrix interface
   ///@{

   //! Returns the total number of non zeros in the matrix.
   virtual int NumNonZeroElems() const;
   /// Gets the columns indexes and values for row *row*.
   /** The return value is always 0 since @a cols and @a srow are copies of the
       values in the matrix. */
   virtual int GetRow(const int row, Array<int> &cols, Vector &srow) const;
   /** @brief If the matrix is square, this method will place 1 on the diagonal
       (i,i) if row i has "almost" zero l1-norm.

       If entry (i,i) does not belong to the sparsity pattern of A, then a error
       will occur. */
   virtual void EliminateZeroRows(const double threshold = 1e-12);

   /// Matrix-Vector Multiplication y = A*x
   virtual void Mult(const Vector & x, Vector & y) const;
   /// Matrix-Vector Multiplication y = y + val*A*x
   virtual void AddMult(const Vector & x, Vector & y, const double val = 1.) const;
   /// MatrixTranspose-Vector Multiplication y = A'*x
   virtual void MultTranspose(const Vector & x, Vector & y) const;
   /// MatrixTranspose-Vector Multiplication y = y + val*A'*x
   virtual void AddMultTranspose(const Vector & x, Vector & y,
                                 const double val = 1.) const;
   ///@}

   /// partial Matrix-Vector Multiplication y = A*x
   void PartMult(const Array<int> &rows, const Vector &x, Vector &y) const;
   /// partial Matrix-Vector Multiplication y = y + a*A*x
   void PartAddMult(const Array<int> &rows, const Vector &x, Vector &y,
                    const double a=1.0) const;

   //! Destructor
   virtual ~BlockMatrix();
   //! If owns_blocks the SparseMatrix objects Aij will be deallocated.
   int owns_blocks;

private:
   //! Given a global row iglobal finds to which row iloc in block iblock belongs to.
   inline void findGlobalRow(int iglobal, int & iblock, int & iloc) const;
   //! Given a global column jglobal finds to which column jloc in block jblock belongs to.
   inline void findGlobalCol(int jglobal, int & jblock, int & jloc) const;

   //! Number of row blocks
   int nRowBlocks;
   //! Number of columns blocks
   int nColBlocks;
   //! row offsets for each block start (length nRowBlocks+1).
   Array<int> row_offsets;
   //! column offsets for each block start (length nColBlocks+1).
   Array<int> col_offsets;
   //! 2D array that stores each block of the BlockMatrix. Aij(iblock, jblock)
   //! == NULL if block (iblock, jblock) is all zeros.
   Array2D<SparseMatrix *> Aij;
};

//! Transpose a BlockMatrix: result = A'
BlockMatrix * Transpose(const BlockMatrix & A);
//! Multiply BlockMatrix matrices: result = A*B
BlockMatrix * Mult(const BlockMatrix & A, const BlockMatrix & B);

inline void BlockMatrix::findGlobalRow(int iglobal, int & iblock,
                                       int & iloc) const
{
   if (iglobal > row_offsets[nRowBlocks])
   {
      mfem_error("BlockMatrix::findGlobalRow");
   }

   for (iblock = 0; iblock < nRowBlocks; ++iblock)
   {
      if (row_offsets[iblock+1] > iglobal) { break; }
   }

   iloc = iglobal - row_offsets[iblock];
}

inline void BlockMatrix::findGlobalCol(int jglobal, int & jblock,
                                       int & jloc) const
{
   if (jglobal > col_offsets[nColBlocks])
   {
      mfem_error("BlockMatrix::findGlobalCol");
   }

   for (jblock = 0; jblock < nColBlocks; ++jblock)
   {
      if (col_offsets[jblock+1] > jglobal) { break; }
   }

   jloc = jglobal - col_offsets[jblock];
}


//! @class SparseBlockMatrix
/**
 * \brief A class storing uniform blocks of given size in sparse fashion
 */
class SparseBlockMatrix : public Matrix
{
public:
   //! Constructor for square sparse block matrices
   /**
     @param s     number of blocks rows and columns
     @param bs    size of the blocks
    */
   SparseBlockMatrix(int s, int bs) : SparseBlockMatrix(s, s, bs, bs) { }

   //! Constructor for rectangular sparse block matrices
   /**
     @param h     number of block rows
     @param w     number of block columns
     @param bh    height of the blocks
     @param bw    width of the blocks
    */
   SparseBlockMatrix(int h, int w, int bh, int bw);

   /// Returns the number of block rows
   inline int NumRowBlocks() const { return blocks.NumRows(); }
   /// Return the number of block columns
   inline int NumColBlocks() const { return blocks.NumCols(); }

   /// Adds the value @a a to the matrix: A(i,j) += a
   /// @note A new block is allocated when being zero
   void Add(const int i, const int j, const double a);

   /// Sets the value @a a to the matrix: A(i,j) = a
   /// @note A new block is allocated when being zero
   void Set(const int i, const int j, const double a);

   /// Return the total number of non-zero blocks
   int NumNonZeroBlocks() const;

   //! Check if block (i,j) is a zero block.
   int IsZeroBlock(int i, int j) const { return (blocks(i,j) == NULL)?(1):(0); }

   /// Returns reference to block (i,j). May be invalid when the block is zero
   DenseMatrix& GetBlock(int i, int j);

   /// Returns reference to constant block (i,j). May be invalid when the block is zero
   const DenseMatrix& GetBlock(int i, int j) const;

   /// @name Matrix interface
   ///@{

   /// Returns reference to a_{ij}.
   virtual double& Elem(int i, int j);

   /// Returns constant reference to a_{ij}.
   virtual const double& Elem(int i, int j) const;

   /// Returns a pointer to (approximation) of the matrix inverse.
   virtual MatrixInverse* Inverse() const
   {
      mfem_error("SparseBlockMatrix::Inverse not implemented \n");
      return static_cast<MatrixInverse*>(NULL);
   }

   /// Prints matrix to stream out.
   virtual void Print(std::ostream & out = mfem::out, int width_ = 4) const;
   ///@}

   ///@name Operator interface
   ///@{

   /// Matrix-Vector Multiplication y = A*x
   virtual void Mult(const Vector &x, Vector &y) const;

   /// MatrixTranspose-Vector Multiplication y = A'*x
   virtual void MultTranspose(const Vector &x, Vector &y) const;
   ///@}

   //! Destructor
   ~SparseBlockMatrix();

private:
   //! 2D array that stores each block of the BlockMatrix. Aij(iblock, jblock)
   //! == NULL if block (iblock, jblock) is all zeros.
   Array2D<DenseMatrix*> blocks;
   //! height of a single block
   int block_height;
   //! width of a single block
   int block_width;

   //! Converts the global to block and sub-block indices
   inline void GetBlockIndices(const int i, const int j, int &ci, int &cj, int &bi,
                               int &bj) const;
};

void SparseBlockMatrix::GetBlockIndices(const int i, const int j, int &ci,
                                        int &cj, int &bi, int &bj) const
{
   ci = i / block_height;
   cj = j / block_width;
   bi = i % block_height;
   bj = j % block_width;
}

}

#endif /* MFEM_BLOCKMATRIX */
