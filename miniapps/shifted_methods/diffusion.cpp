﻿//                       MFEM
//
// Compile with: make diffusion
//
// Sample runs:
// mpirun -np 1 diffusion -m ../../data/inline-quad.mesh  -rs 0 -vis -o 2
// mpirun -np 1 diffusion -m quad.mesh -rs 0 -o 2 -st 1 -ex -lst 1
// mpirun -np 1 diffusion -m quad.mesh -rs 0 -o 2 -st 1 -ex -lst 3
// mpirun -np 1 diffusion -m quad.mesh -rs 0 -o 2 -st 1 -ex -lst 4
// mpirun -np 4 diffusion -m quad.mesh -rs 2 -o 3 -st 1 -no-ex -lst 5
#include "../../mfem.hpp"
#include <fstream>
#include <iostream>
#include "sbm-aux.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   /// 1. Parse command-line options.
   const char *mesh_file = "../data/square-disc.mesh";
   int order = 2;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   int ser_ref_levels = 0;
   bool exact = true;
   int level_set_type = 1;
   int ho_terms = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&exact, "-ex", "--exact", "-no-ex",
                  "--no-exact",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&level_set_type, "-lst", "--level-set-type",
                  "level-set-type:");
   args.AddOption(&ho_terms, "-ho", "--high-order",
                  "Additional high-order terms to include");
   /// IMPORTANT: Level set type notes
   /// 1 - circlular hole of radius 0.2 at the center of domain [0, 1].
   ///     -nabla^2 = 1. Exact solution can be generated by using a body-fitted
   ///     fine mesh.
   ///
   /// 2 - circular hole of radius 0.2 at the center of domain [0, 1]. Solution
   /// is of the form y=x^p+y^p. We use this to make sure we get exact solution.
   /// p is currently set to 2. It can be changed in sbm-aux.hpp.
   ///
   /// 3. - Walls are at y = 0 to 1. In this case we stretch the mesh from
   /// [0,1] to [-1e-4, 1] to get shifted faces. The solution is analytic
   /// sinusoidal function. - Used for convergence study.

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();
   for (int lev = 0; lev < ser_ref_levels; lev++) { mesh.UniformRefinement(); }

   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
   {
      int par_ref_levels = 0;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh.UniformRefinement();
      }
   }

   // 4. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
   }
   ParFiniteElementSpace pfespace(&pmesh, fec);

   Vector vxyz;

   ParFiniteElementSpace pfespace_mesh(&pmesh, fec, dim);
   pmesh.SetNodalFESpace(&pfespace_mesh);
   ParGridFunction x_mesh(&pfespace_mesh);
   pmesh.SetNodalGridFunction(&x_mesh);
   vxyz = *pmesh.GetNodes();
   int nodes_cnt = vxyz.Size()/dim;
   if (level_set_type == 3) { //stretch quadmesh from [0, 1] to [-1.e-4, 1]
       for (int i = 0; i < nodes_cnt; i++) {
           vxyz(i+nodes_cnt) = (1.+1.e-4)*vxyz(i+nodes_cnt)-1.e-4;
       }
   }
   pmesh.SetNodes(vxyz);
   pfespace.ExchangeFaceNbrData();

   // Setup FESpace for L2 function (used to mark element flags etc.)
   L2_FECollection fecl2 = L2_FECollection(0, dim);
   L2_FECollection fecl2ho = L2_FECollection(order, dim);
   ParFiniteElementSpace pfesl2(&pmesh, &fecl2);
   ParFiniteElementSpace pfesl2ho(&pmesh, &fecl2ho);
   cout << "Number of finite element unknowns: "
        << pfespace.GetTrueVSize() << endl;

   // 5. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction x(&pfespace);
   // ParGridFunction for level_set_value.
   ParGridFunction level_set_val(&pfespace);

   // 6. Determine if each element in the ParMesh is inside the actual domain,
   //    partially cut by the actual domain boundary, or completely outside
   //    the domain.
   Dist_Level_Set_Coefficient dist_fun_level_coef(level_set_type);
   level_set_val.ProjectCoefficient(dist_fun_level_coef);
   level_set_val.ExchangeFaceNbrData();

   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   FaceElementTransformations *tr = NULL;

   // Set trim flag based on the distance field
   // 0 if completely in the domain
   // 1 if completely outside the domain
   // 2 if partially inside the domain
   ParGridFunction elem_marker(&pfesl2);
   Array<int> trim_flag(pmesh.GetNE()+pmesh.GetNSharedFaces());
   trim_flag = 0;
   Vector vals;
   // Check elements on the current MPI rank
   for (int i = 0; i < pmesh.GetNE(); i++)
   {
      ElementTransformation *Tr = pmesh.GetElementTransformation(i);
      const IntegrationRule &ir =
         IntRulesLo.Get(pmesh.GetElementBaseGeometry(i), 4*Tr->OrderJ());
      level_set_val.GetValues(i, ir, vals);

      int count = 0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         double val = vals(j);
         if (val <= 0.) { count++; }
      }
      if (count == ir.GetNPoints()) // completely outside
      {
         trim_flag[i] = 1;
      }
      else if (count > 0) // partially outside
      {
         trim_flag[i] = 2;
      }
   }

   // Check neighbors on the adjacent MPI rank
   for (int i = pmesh.GetNE(); i < pmesh.GetNE()+pmesh.GetNSharedFaces(); i++)
   {
       int shared_fnum = i-pmesh.GetNE();
       tr = pmesh.GetSharedFaceTransformations(shared_fnum);
       int Elem2NbrNo = tr->Elem2No - pmesh.GetNE();

       ElementTransformation *eltr =
               pfespace.GetFaceNbrElementTransformation(Elem2NbrNo);
       const IntegrationRule &ir =
         IntRulesLo.Get(pfespace.GetFaceNbrFE(Elem2NbrNo)->GetGeomType(),
                        4*eltr->OrderJ());

       const int nip = ir.GetNPoints();
       vals.SetSize(nip);
       int count = 0;
       for (int j = 0; j < nip; j++) {
          const IntegrationPoint &ip = ir.IntPoint(j);
          vals[j] = level_set_val.GetValue(tr->Elem2No, ip);
          if (vals[j] <= 0.) { count++; }
       }

      if (count == ir.GetNPoints()) // completely outside
      {
         trim_flag[i] = 1;
      }
      else if (count > 0) // partially outside
      {
         trim_flag[i] = 2;
      }
   }

   for (int i = 0; i < elem_marker.Size(); i++)
   {
      elem_marker(i) = trim_flag[i]*1.;
   }

   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock << "solution\n" << pmesh << elem_marker << flush;
      sol_sock << "window_title 'Element flags'\n"
               << "window_geometry "
               << 0 << " " << 0 << " " << 350 << " " << 350 << "\n"
               << "keys Rjmpc" << endl;
   }

   // 7. Determine Shifted boundary faces of the mesh and get a list of dofs
   //    associated with those faces.
   Array<int> sbm_dofs; // Array of dofs on sbm faces
   Array<int> dofs; // work array

   // First we check interior faces of the mesh (excluding interior faces that
   // are on the processor boundaries)
   double count1 = 0;
   for (int i = 0; i < pmesh.GetNumFaces(); i++)
   {
      FaceElementTransformations *tr = NULL;
      tr = pmesh.GetInteriorFaceTransformations (i);
      const int faceno = i;
      if (tr != NULL)
      {
         count1 += 1;
         int ne1 = tr->Elem1No;
         int ne2 = tr->Elem2No;
         int te1 = trim_flag[ne1], te2 = trim_flag[ne2];
         if (te1 == 2 && te2 == 0)
         {
             pfespace.GetFaceDofs(faceno, dofs);
             sbm_dofs.Append(dofs);
         }
         if (te1 == 0 && te2 == 2)
         {
             pfespace.GetFaceDofs(faceno, dofs);
             sbm_dofs.Append(dofs);
         }
      }
   }


   // Here we add boundary faces that we want to model as SBM faces.
   // For the method where we clip inside the domain, a boundary face
   // has to be set as SBM face using its attribute.
   double count2 = 0;
   for (int i = 0; i < pmesh.GetNBE(); i++)
   {
      int attr = pmesh.GetBdrAttribute(i);
      FaceElementTransformations *tr;
      tr = pmesh.GetBdrFaceTransformations (i);
      if (tr != NULL) {
          if (attr == 100) { // add all boundary faces with attr=100 as SBM faces
              count2 += 1;
              int ne1 = tr->Elem1No;
              int te1 = trim_flag[ne1];
              const int faceno = pmesh.GetBdrFace(i);
              if (te1 == 0)
              {
                 pfespace.GetFaceDofs(faceno, dofs);
                 sbm_dofs.Append(dofs);
              }
          }
      }
   }

   // Now we add interior faces that are on processor boundaries.
   double count3 = 0;
   double count3b = 0;
   for (int i = 0; i < pmesh.GetNSharedFaces(); i++)
   {
      tr = pmesh.GetSharedFaceTransformations(i);
      if (tr != NULL)
      {
         count3b += 1;
         int ne1 = tr->Elem1No;
         int te1 = trim_flag[ne1];
         int te2 = trim_flag[i+pmesh.GetNE()];
         const int faceno = pmesh.GetSharedFace(i);
         // Add if the element on this proc is completely inside the domain
         // and the the element on other proc is not
         if (te2 == 2 && te1 == 0)
         {
            count3 += 1;
            pfespace.GetFaceDofs(faceno, dofs);
            sbm_dofs.Append(dofs);
         }
      }
   }


   // 8. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    To do this, we first make a list of all dofs that are on the real boundary
   //    of the mesh, then add all the dofs of the elements that are completely
   //    outside or intersect shifted boundary. Then we remove the dofs from
   //    SBM faces.

   // Make a list of dofs on all boundaries
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(pmesh.bdr_attributes.Max());
   if (pmesh.bdr_attributes.Size())
   {
      ess_bdr = 1;
   }
   Array<int> ess_vdofs_bdr;
   pfespace.GetEssentialVDofs(ess_bdr, ess_vdofs_bdr);

   // Get all dofs associated with elements outside the domain or intersected
   // by the boundary.
   Array<int> ess_vdofs_hole(ess_vdofs_bdr.Size());
   ess_vdofs_hole = 0;
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      if (trim_flag[e] > 0)
      {
         pfespace.GetElementVDofs(e, dofs);
         for (int i = 0; i < dofs.Size(); i++) {
             ess_vdofs_hole[dofs[i]] = -1;
         }
      }
   }

   // Combine the lists to mark essential dofs.
   for (int i = 0; i < ess_vdofs_hole.Size(); i++)
   {
      if (ess_vdofs_bdr[i] == -1) { ess_vdofs_hole[i] = -1; }
   }

   // Unmark dofs that are on SBM faces (but not on dirichlet boundaries)
   for (int i = 0; i < sbm_dofs.Size(); i++) {
       if (ess_vdofs_bdr[sbm_dofs[i]] != -1) { 
          ess_vdofs_hole[sbm_dofs[i]] = 0;
       }
   }

   // Synchronize
   for (int i = 0; i < ess_vdofs_hole.Size() ; i++) {
       ess_vdofs_hole[i] += 1;
   }

   pfespace.Synchronize(ess_vdofs_hole);

   for (int i = 0; i < ess_vdofs_hole.Size() ; i++) {
       ess_vdofs_hole[i] -= 1;
   }

   // convert to tdofs
   Array<int> ess_tdofs;
   pfespace.GetRestrictionMatrix()->BooleanMult(ess_vdofs_hole,
                                                ess_tdofs);
   pfespace.MarkerToList(ess_tdofs, ess_tdof_list);


   // Compute Distance Vector - Use analytic distance vectors for now.
   auto distance_vec_space = new ParFiniteElementSpace(pfespace.GetParMesh(),
                                                       pfespace.FEColl(), dim);
   ParGridFunction distance(distance_vec_space);

   // 9. Get the Distance from the level set either using a numerical approach
   //    or project an exact analytic function.
//   HeatDistanceSolver dist_func(1.0);
//   dist_func.smooth_steps = 1;
//   dist_func.ComputeVectorDistance(dist_fun_level_coef, distance);

   VectorCoefficient *dist_vec = NULL;
   if (true) {
       Dist_Vector_Coefficient *dist_vec_fcoeff =
               new Dist_Vector_Coefficient(dim, level_set_type);
       dist_vec = dist_vec_fcoeff;
       distance.ProjectDiscCoefficient(*dist_vec);
   }

   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock << "solution\n" << pmesh << distance << flush;
      sol_sock << "window_title 'Distance Vector'\n"
               << "window_geometry "
               << 350 << " " << 350 << " " << 350 << " " << 350 << "\n"
               << "keys Rjmpcvv" << endl;
   }

   // 12. SBM integration parameter - alpha
   double alpha = 1;

   // 13. Set up the linear form b(.) which corresponds to the right-hand side of
   //     the FEM linear system.
   ParLinearForm b(&pfespace);
   SBMFunctionCoefficient rhs_f(rhs_fun, level_set_type);
   b.AddDomainIntegrator(new DomainLFIntegrator(rhs_f), trim_flag);

   SBMFunctionCoefficient dbcCoef(dirichlet_velocity, level_set_type);
   b.AddShiftedBdrFaceIntegrator(new SBM2LFIntegrator(dbcCoef, alpha, *dist_vec, ho_terms),
                                 trim_flag);
   b.Assemble();

   // 14. Set up the bilinear form a(.,.) on the finite element space
   //     corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //     domain integrator and SBM integrator.
   ParBilinearForm a(&pfespace);
   ConstantCoefficient one(1.);

   a.AddDomainIntegrator(new DiffusionIntegrator(one), trim_flag);
   a.AddShiftedBdrFaceIntegrator(new SBM2Integrator(alpha, *dist_vec, ho_terms),
                                 trim_flag);

   // 15. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations.
   a.Assemble();

   // Project the exact solution as an initial condition for dirichlet boundaries.
   x = 0;
   x.ProjectCoefficient(dbcCoef);
   // Zero out non-essential boundaries.
   for (int i = 0; i < ess_vdofs_hole.Size(); i++) {
       if (ess_vdofs_hole[i] != -1) { x(i) = 0.; }
   }

   // 16. Form the linear system and solve it.
   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   Solver *S = NULL;
   Solver *prec = NULL;
   prec = new HypreBoomerAMG;
   BiCGSTABSolver *bicg = new BiCGSTABSolver(MPI_COMM_WORLD);
   bicg->SetRelTol(1e-12);
   bicg->SetMaxIter(2000);
   bicg->SetPrintLevel(1);
   bicg->SetPreconditioner(*prec);
   bicg->SetOperator(*A);
   S = bicg;
   S->Mult(B, X);

   // 17. Recover the solution as a finite element grid function.
   a.RecoverFEMSolution(X, b, x);

   // 18. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("ex1-sbm.mesh");
   mesh_ofs.precision(8);
   pmesh.PrintAsOne(mesh_ofs);
   ofstream sol_ofs("ex1-sbm.gf");
   sol_ofs.precision(8);
   x.SaveAsOne(sol_ofs);

   // 19. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock << "solution\n" << pmesh << x << flush;
      sol_sock << "window_title 'Solution'\n"
               << "window_geometry "
               << 350 << " " << 0 << " " << 350 << " " << 350 << "\n"
               << "keys Rj" << endl;
   }

   // 20. Construct an error gridfunction if the exact solution is known.
   ParGridFunction err(x);
   Vector pxyz(dim);
   pxyz(0) = 0.;
   for (int i = 0; i < nodes_cnt; i++) {
       pxyz(0) = vxyz(i);;
       pxyz(1) = vxyz(i+nodes_cnt);
       double exact_val = dirichlet_velocity(pxyz, level_set_type);
       err(i) = std::fabs(x(i) - exact_val);
   }

   double global_error = err.Norml2();
   if (myid == 0 && level_set_type != 1) {
       std::cout << global_error << " Global error - L2 norm.\n";
   }

   if (visualization && level_set_type >= 3 && level_set_type <= 4)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock << "solution\n" << pmesh << err << flush;
      sol_sock << "window_title 'Error'\n"
               << "window_geometry "
               << 700 << " " << 0 << " " << 350 << " " << 350 << "\n"
               << "keys Rj" << endl;
   }

   // 15. Free the used memory.
   delete prec;
   delete S;
   delete dist_vec;
   delete distance_vec_space;
   delete fec;

   MPI_Finalize();

   return 0;
}
