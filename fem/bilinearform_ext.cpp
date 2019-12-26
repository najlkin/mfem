// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

// Implementations of classes FABilinearFormExtension, EABilinearFormExtension,
// PABilinearFormExtension and MFBilinearFormExtension.

#include "../general/forall.hpp"
#include "bilinearform.hpp"

namespace mfem
{

BilinearFormExtension::BilinearFormExtension(BilinearForm *form)
   : Operator(form->Size()), a(form)
{
   // empty
}

const Operator *BilinearFormExtension::GetProlongation() const
{
   return a->GetProlongation();
}

const Operator *BilinearFormExtension::GetRestriction() const
{
   return a->GetRestriction();
}


// Data and methods for element-assembled bilinear forms
EABilinearFormExtension::EABilinearFormExtension(BilinearForm *form)
   : BilinearFormExtension(form),
     fes(a->FESpace())
{
   elem_restrict_nat = fes->GetElementRestriction(
                          ElementDofOrdering::NATIVE);
   if (elem_restrict_nat)
   {
      localX.SetSize(elem_restrict_nat->Height(), Device::GetMemoryType());
      localY.SetSize(elem_restrict_nat->Height(), Device::GetMemoryType());
      localY.UseDevice(true); // ensure 'localY = 0.0' is done on device
   }
}

void EABilinearFormExtension::Assemble()
{
   a->ComputeElementMatrices();
}

void EABilinearFormExtension::Update()
{
   fes = a->FESpace();
   height = width = fes->GetVSize();
   elem_restrict_nat = fes->GetElementRestriction(
                          ElementDofOrdering::NATIVE);
   if (elem_restrict_nat)
   {
      localX.SetSize(elem_restrict_nat->Height());
      localY.SetSize(elem_restrict_nat->Height());
   }
}

void EABilinearFormExtension::FormSystemMatrix(const Array<int> &ess_tdof_list,
                                               OperatorHandle &A)
{
   const Operator* P = fes->GetProlongationMatrix();
   Operator *rap = this;
   if (P) { rap = new RAPOperator(*P, *this, *P); }
   const bool own_A = (rap!=this);
   A.Reset(new ConstrainedOperator(rap, ess_tdof_list, own_A));
}

void EABilinearFormExtension::FormLinearSystem(const Array<int> &ess_tdof_list,
                                               Vector &x, Vector &b,
                                               OperatorHandle &A,
                                               Vector &X, Vector &B,
                                               int copy_interior)
{
   Operator *oper;
   Operator::FormLinearSystem(ess_tdof_list, x, b, oper, X, B, copy_interior);
   A.Reset(oper); // A will own oper
}

void EABilinearFormExtension::Mult(const Vector &x, Vector &y) const
{
   const int NE = fes->GetNE();
   DenseTensor tensorX, tensorY;
   const int ndofs_el = fes->GetFE(0)->GetDof();
   const int vdim = fes->GetVDim();

   if (elem_restrict_nat)
   {
      MFEM_ASSERT(elem_restrict_nat->Height() == ndofs_el * vdim * NE,
                  "E-vector size is not NDOFS x VDIM x NE");
      tensorX.UseExternalData(localX.GetData(), ndofs_el, vdim, NE);
      tensorY.UseExternalData(localY.GetData(), ndofs_el, vdim, NE);

      elem_restrict_nat->Mult(x, localX);
   }
   else
   {
      MFEM_ASSERT(x.Size() == ndofs_el * vdim * NE,
                  "E-vector size is not NDOFS x VDIM x NE");
      tensorX.UseExternalData(x.GetData(), ndofs_el, vdim, NE);
      tensorY.UseExternalData(y.GetData(), ndofs_el, vdim, NE);
   }

   DenseMatrix elmat;

   for (int i = 0; i < NE; i++)
   {
      a->ComputeElementMatrix(i, elmat);

      mfem::Mult(elmat, tensorX(i), tensorY(i));
   }

   if (elem_restrict_nat)
   {
      elem_restrict_nat->MultTranspose(localY, y);
   }
}

void EABilinearFormExtension::MultTranspose(const Vector &x, Vector &y) const
{
   const int NE = fes->GetNE();
   DenseTensor tensorX, tensorY;
   const int ndofs_el = fes->GetFE(0)->GetDof();
   const int vdim = fes->GetVDim();

   if (elem_restrict_nat)
   {
      MFEM_ASSERT(elem_restrict_nat->Height() == ndofs_el * vdim * NE,
                  "E-vector size is not NDOFS x VDIM x NE");
      tensorX.UseExternalData(localX.GetData(), ndofs_el, vdim, NE);
      tensorY.UseExternalData(localY.GetData(), ndofs_el, vdim, NE);

      elem_restrict_nat->Mult(x, localX);
   }
   else
   {
      MFEM_ASSERT(x.Size() == ndofs_el * vdim * NE,
                  "E-vector size is not NDOFS x VDIM x NE");
      tensorX.UseExternalData(x.GetData(), ndofs_el, vdim, NE);
      tensorY.UseExternalData(y.GetData(), ndofs_el, vdim, NE);
   }

   DenseMatrix elmat;

   for (int i = 0; i < NE; i++)
   {
      a->ComputeElementMatrix(i, elmat);

      mfem::MultTranspose(elmat, tensorX(i), tensorY(i));
   }

   if (elem_restrict_nat)
   {
      elem_restrict_nat->MultTranspose(localY, y);
   }
}

// Data and methods for partially-assembled bilinear forms
PABilinearFormExtension::PABilinearFormExtension(BilinearForm *form)
   : BilinearFormExtension(form),
     trialFes(a->FESpace()), testFes(a->FESpace())
{
   elem_restrict_lex = trialFes->GetElementRestriction(
                          ElementDofOrdering::LEXICOGRAPHIC);
   if (elem_restrict_lex)
   {
      localX.SetSize(elem_restrict_lex->Height(), Device::GetMemoryType());
      localY.SetSize(elem_restrict_lex->Height(), Device::GetMemoryType());
      localY.UseDevice(true); // ensure 'localY = 0.0' is done on device
   }
}

void PABilinearFormExtension::Assemble()
{
   Array<BilinearFormIntegrator*> &integrators = *a->GetDBFI();
   const int integratorCount = integrators.Size();
   for (int i = 0; i < integratorCount; ++i)
   {
      integrators[i]->AssemblePA(*a->FESpace());
   }
}

void PABilinearFormExtension::Update()
{
   FiniteElementSpace *fes = a->FESpace();
   height = width = fes->GetVSize();
   trialFes = fes;
   testFes = fes;
   elem_restrict_lex = trialFes->GetElementRestriction(
                          ElementDofOrdering::LEXICOGRAPHIC);
   if (elem_restrict_lex)
   {
      localX.SetSize(elem_restrict_lex->Height());
      localY.SetSize(elem_restrict_lex->Height());
   }
}

void PABilinearFormExtension::FormSystemMatrix(const Array<int> &ess_tdof_list,
                                               OperatorHandle &A)
{
   const Operator* trialP = trialFes->GetProlongationMatrix();
   const Operator* testP  = testFes->GetProlongationMatrix();
   Operator *rap = this;
   if (trialP) { rap = new RAPOperator(*testP, *this, *trialP); }
   const bool own_A = (rap!=this);
   A.Reset(new ConstrainedOperator(rap, ess_tdof_list, own_A));
}

void PABilinearFormExtension::FormLinearSystem(const Array<int> &ess_tdof_list,
                                               Vector &x, Vector &b,
                                               OperatorHandle &A,
                                               Vector &X, Vector &B,
                                               int copy_interior)
{
   Operator *oper;
   Operator::FormLinearSystem(ess_tdof_list, x, b, oper, X, B, copy_interior);
   A.Reset(oper); // A will own oper
}

void PABilinearFormExtension::Mult(const Vector &x, Vector &y) const
{
   Array<BilinearFormIntegrator*> &integrators = *a->GetDBFI();

   const int iSz = integrators.Size();
   if (elem_restrict_lex)
   {
      elem_restrict_lex->Mult(x, localX);
      localY = 0.0;
      for (int i = 0; i < iSz; ++i)
      {
         integrators[i]->AddMultPA(localX, localY);
      }
      elem_restrict_lex->MultTranspose(localY, y);
   }
   else
   {
      y.UseDevice(true); // typically this is a large vector, so store on device
      y = 0.0;
      for (int i = 0; i < iSz; ++i)
      {
         integrators[i]->AddMultPA(x, y);
      }
   }
}

void PABilinearFormExtension::MultTranspose(const Vector &x, Vector &y) const
{
   Array<BilinearFormIntegrator*> &integrators = *a->GetDBFI();
   const int iSz = integrators.Size();
   if (elem_restrict_lex)
   {
      elem_restrict_lex->Mult(x, localX);
      localY = 0.0;
      for (int i = 0; i < iSz; ++i)
      {
         integrators[i]->AddMultTransposePA(localX, localY);
      }
      elem_restrict_lex->MultTranspose(localY, y);
   }
   else
   {
      y.UseDevice(true);
      y = 0.0;
      for (int i = 0; i < iSz; ++i)
      {
         integrators[i]->AddMultTransposePA(x, y);
      }
   }
}

} // namespace mfem
