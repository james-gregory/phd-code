// This header file
// contains a bunch of derived classes which are easier to find when separated,
// we have derived elements whcih output some specified quantities, a transvsersely
// isotropic strain energy function constitutive law, several strain energy functions
// and the definition of a displaced boundary

//====================================================================
//Continuous Pressure QPVDElements
//====================================================================

namespace oomph
{
 template<unsigned DIM>
 class MyQPVDElementWithContinuousPressure : public virtual QPVDElementWithContinuousPressure<DIM>
{
public:

 /// Constructor: Call constructor of underlying element
  MyQPVDElementWithContinuousPressure() :  QPVDElementWithContinuousPressure<DIM>() {}

  /// Overload output function
  void output(std::ostream &outfile, const unsigned &n_plot)
   {
    Vector<double> s(3);
    Vector<double> x(3);
    Vector<double> xi(3);
    DenseMatrix<double> sigma(3,3);
    DenseMatrix<double> deformation_gradient(3,3);
    DenseMatrix<double> cauchy_stress(3,3);
    DenseMatrix<double> right_cauchy_green(3,3);
    Vector<double> M(3);
    Vector<double> FM(3);

   // initialise to zero
   for(unsigned i = 0; i < 3; i++)
   {
 	  M[i] = 0;
 	  FM[i] = 0;
   }

    //Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << ", K=" << n_plot << std::endl;

      //Loop over element nodes
      for(unsigned l3=0;l3<n_plot;l3++)
       {
        s[2] = -1.0 + l3*2.0/(n_plot-1);
        for(unsigned l2=0;l2<n_plot;l2++)
         {
          s[1] = -1.0 + l2*2.0/(n_plot-1);
          for(unsigned l1=0;l1<n_plot;l1++)
           {
            s[0] = -1.0 + l1*2.0/(n_plot-1);

        // Get Eulerian and Lagrangian coordinates and the stress
        this->interpolated_x(s,x);
        this->interpolated_xi(s,xi);
        this->get_stress(s,sigma);

        //Output the x,y,z coordinates
        for(unsigned i=0;i<3;i++)
   	 { outfile << x[i] << " ";}

 	 // number of local nodes for the element
 	 unsigned n_node = this->nnode();
 	 // number of coordinate types required to interpolate between the nodes
 	 unsigned n_position_type = this->nnodal_position_type();
 	 // create the shape functions
 	 Shape psi(n_node,n_position_type);
 	 // and their derivatives
 	 DShape dpsi(n_node,n_position_type,3);
 	 // calculate the shape functions
 	 this->shape(s,psi);
 	 // and their derivatives
 	 this->dshape_lagrangian(s,psi,dpsi);

 	 // initialise to zero
 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 deformation_gradient(i,j) = 0.0;
 			 cauchy_stress(i,j) = 0.0;
 			 right_cauchy_green(i,j) = 0.0;
 		 }
 	 }

 	 // loop over the number of nodes
 	 for(unsigned j=0;j<n_node;j++)
 	 {
 		 // loop over the dimensions
 		 for(unsigned i=0;i<3;i++)
 		 {
 			 // loop over the derivative directions
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 deformation_gradient(i,k) += this->nodal_position(j,i)*dpsi(j,k);
 			 }
 		 }
 	 }

 	 // calculate the cauchy stress from the second piola kirchhoff stress (sigma)
 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 for(unsigned l=0;l<3;l++)
 				 {
 					 cauchy_stress(i,j) +=
 					 deformation_gradient(i,k)*sigma(k,l)*deformation_gradient(j,l);
 				 }
 			 }
 		 }
 	 }


 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 right_cauchy_green(i,j) += deformation_gradient(k,i)*deformation_gradient(k,j);
 			 }
 		 }
 	 }

         //Output xi[0], xi[1], xi[2]
         for(unsigned i=0;i<3;i++)
         {outfile << xi[i] << " ";}

		 M[0] = 0.0;
		 M[1] = 0.0;
		 M[2] = 1.0;

 		// compute the deformed fibre direction
 		for(unsigned i = 0; i < 3; i++)
 		{
 			for(unsigned j = 0; j < 3; j++)
 			{
 				FM[i] = deformation_gradient(i,j)*M[j];
 			}
 		}
 		// turn into unit vector
 		double FM_magnitude = sqrt(FM[0]*FM[0] + FM[1]*FM[1] + FM[2]*FM[2]);
 		for(unsigned i = 0; i < 3; i++)
 		{
 			FM[i] = FM[i]/FM_magnitude;
 		}


 		// calculate I4
 		double I4 = 0.0;
 		for(unsigned i=0;i<3;i++)
 		  {
 			for(unsigned j=0;j<3;j++)
 			  {
 				I4 += M[i]*right_cauchy_green(i,j)*M[j];
 			  }
 		  }

        //Output stress
        outfile    << sigma(0,0) << " "
        	       << sigma(1,1) << " "
        	       << sigma(2,2) << " "
        	       << sigma(0,1) << " "
        	       << sigma(0,2) << " "
        	       << sigma(1,2) << std::endl;


	}
}
}
} // end of output



};

//=========================================================================
///FaceGeometry of the 3D QPVDElementWithContinuousPressure
//=========================================================================
template<>
class FaceGeometry<MyQPVDElementWithContinuousPressure<3> >:
 public virtual SolidQElement<2,3>
{
  public:
 //Make sure that we call the constructor of the SolidQElement
 //Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2,3>() {}
};

//===============================================================
/// FaceGeometry of FaceGeometry
/// for 3D MyQPVDElementWithContinuousPressure element
//===============================================================
template<>
class FaceGeometry<FaceGeometry<MyQPVDElementWithContinuousPressure<3> > >:
public virtual SolidQElement<1,3>
{
  public:
 /// Constructor must call constructor of the underlying element
  FaceGeometry() : SolidQElement<1,3>() {}
};

//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

  template<unsigned DIM>
  class MyRefineableQPVDElementWithContinuousPressure : public virtual RefineableQPVDElementWithContinuousPressure<DIM>
{

public:

 /// Constructor: Call constructor of underlying element
  MyRefineableQPVDElementWithContinuousPressure() :  RefineableQPVDElementWithContinuousPressure<DIM>() {}

  /// Overload output function
  void output(std::ostream &outfile, const unsigned &n_plot)
   {
	Vector<double> s(3);
	Vector<double> x(3);
	Vector<double> xi(3);
	DenseMatrix<double> sigma(3,3);
	DenseMatrix<double> deformation_gradient(3,3);
	DenseMatrix<double> cauchy_stress(3,3);
	DenseMatrix<double> right_cauchy_green(3,3);
	Vector<double> M(3);
	Vector<double> FM(3);



	//Tecplot header info
	outfile << "ZONE I=" << n_plot << ", J=" << n_plot << ", K=" << n_plot << std::endl;

	  //Loop over element nodes
	  for(unsigned l3=0;l3<n_plot;l3++)
	   {
		s[2] = -1.0 + l3*2.0/(n_plot-1);
		for(unsigned l2=0;l2<n_plot;l2++)
		 {
		  s[1] = -1.0 + l2*2.0/(n_plot-1);
		  for(unsigned l1=0;l1<n_plot;l1++)
		   {
			s[0] = -1.0 + l1*2.0/(n_plot-1);

		// Get Eulerian and Lagrangian coordinates and the stress
		this->interpolated_x(s,x);
		this->interpolated_xi(s,xi);
		this->get_stress(s,sigma);

		//Output the x,y,z coordinates
	    for(unsigned i=0;i<3;i++)
	   { outfile << x[i] << " ";}

   // number of local nodes for the element
   unsigned n_node = this->nnode();
   // number of coordinate types required to interpolate between the nodes
   unsigned n_position_type = this->nnodal_position_type();
   // create the shape functions
   Shape psi(n_node,n_position_type);
   // and their derivatives
   DShape dpsi(n_node,n_position_type,3);
   // calculate the shape functions
   this->shape(s,psi);


   // and their derivatives
   this->dshape_lagrangian(s,psi,dpsi);

   // initialise to zero
   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   deformation_gradient(i,j) = 0.0;
		   cauchy_stress(i,j) = 0.0;
		   right_cauchy_green(i,j) = 0.0;
	   }
   }


   // initialise to zero
   for(unsigned i = 0; i < 3; i++)
   {
	M[i] = 0;
	FM[i] = 0;
   }

   // loop over the number of nodes
   for(unsigned j=0;j<n_node;j++)
   {
	   // loop over the dimensions
	   for(unsigned i=0;i<3;i++)
	   {
		   // loop over the derivative directions
		   for(unsigned k=0;k<3;k++)
		   {
			   deformation_gradient(i,k) += this->nodal_position(j,i)*dpsi(j,k);
		   }
	   }
   }

   // calculate the cauchy stress from the second piola kirchhoff stress (sigma)
   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   for(unsigned k=0;k<3;k++)
		   {
			   for(unsigned l=0;l<3;l++)
			   {
				   cauchy_stress(i,j) +=
				   deformation_gradient(i,k)*sigma(k,l)*deformation_gradient(j,l);
			   }
		   }
	   }
   }


   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   for(unsigned k=0;k<3;k++)
		   {
			   right_cauchy_green(i,j) += deformation_gradient(k,i)*deformation_gradient(k,j);
		   }
	   }
   }

		 //Output xi[0], xi[1], xi[2]
		 for(unsigned i=0;i<3;i++)
		 {outfile << xi[i] << " ";}

	   M[0] = 0.0;
	   M[1] = 0.0;
	   M[2] = 1.0;


	  // compute the deformed fibre direction
	  for(unsigned i = 0; i < 3; i++)
	  {
		  for(unsigned j = 0; j < 3; j++)
		  {
			  FM[i] = deformation_gradient(i,j)*M[j];
		  }
	  }
	  // turn into unit vector
	  double FM_magnitude = sqrt(FM[0]*FM[0] + FM[1]*FM[1] + FM[2]*FM[2]);
	  for(unsigned i = 0; i < 3; i++)
	  {
		  FM[i] = FM[i]/FM_magnitude;
	  }


	  // calculate I4
	  double I4 = 0.0;
	  for(unsigned i=0;i<3;i++)
		{
		  for(unsigned j=0;j<3;j++)
			{
			  I4 += M[i]*right_cauchy_green(i,j)*M[j];
			}
		}

		//Output stress
        outfile    << sigma(0,0) << " "
        	       << sigma(1,1) << " "
        	       << sigma(2,2) << " "
        	       << sigma(0,1) << " "
        	       << sigma(0,2) << " "
        	       << sigma(1,2) << std::endl;

	}
}
}
} // end of output



  }; //end MyRefineableQPVDElementWithContinuousPressure

//=========================================================================
///FaceGeometry of the 3D MyRefineableQPVDElementWithContinuousPressure
//=========================================================================
template<>
class FaceGeometry<MyRefineableQPVDElementWithContinuousPressure<3> >:
 public virtual SolidQElement<2,3>
{
  public:
 //Make sure that we call the constructor of the SolidQElement
 //Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2,3>() {}
};

//===============================================================
/// FaceGeometry of FaceGeometry
/// for 3D MyRefineableQPVDElementWithContinuousPressure element
//===============================================================
template<>
class FaceGeometry<FaceGeometry<MyRefineableQPVDElementWithContinuousPressure<3> > >:
public virtual SolidQElement<1,3>
{
  public:
 /// Constructor must call constructor of the underlying element
  FaceGeometry() : SolidQElement<1,3>() {}
};


//===============================================================
// A class derived from ImposeDisplacementByLagrangeMultiplierElements
// that include a new output function designed to output the
// resultant force over the top of the cylinder
//===============================================================
template<class ELEMENT>
class MyLagrangeMultiplierElement : public ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>
{
public:


  MyLagrangeMultiplierElement(FiniteElement* const &element_pt, const int &face_index):
    ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(element_pt,face_index) { }






		/// \short Output function
	  void output(std::ostream &outfile, const unsigned &n_plot)
	   {
	    // Elemental dimension
	    unsigned dim_el=this->dim();

	    //Find the number of positional types
	    unsigned n_position_type = this->nnodal_position_type();

	 #ifdef PARANOID
	    if(n_position_type!=1)
	     {
	      throw OomphLibError(
	       "ImposeDisplacementByLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
	       OOMPH_CURRENT_FUNCTION,
	       OOMPH_EXCEPTION_LOCATION);
	     }
	 #endif


	    //Local coord
	    Vector<double> s(dim_el);

	    // # of nodes,
	    unsigned n_node=this->nnode();
	    Shape psi(n_node,n_position_type);

	    // Tecplot header info
	    //outfile << this->tecplot_zone_string(n_plot);

	    // Loop over plot points
	    unsigned num_plot_points=this->nplot_points(n_plot);
	    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
	     {
	      // Get local coordinates of plot point
	      this->get_s_plot(iplot,n_plot,s);

	      // Get shape function
	      this->shape(s,psi);

	      //Calculate the Eulerian coordinates and Lagrange multiplier
	      Vector<double> x(dim_el+1,0.0);
	      Vector<double> lambda(dim_el+1,0.0);
	      Vector<double> zeta(dim_el,0.0);
//				Vector<double> resultant_force(dim_el+1,0.0);
	      for(unsigned j=0;j<n_node;j++)
	       {
	        // Cast to a boundary node
	        BoundaryNodeBase *bnod_pt =
	         dynamic_cast<BoundaryNodeBase*>(this->node_pt(j));

	        // get the node pt
	        Node* nod_pt = this->node_pt(j);

	        // Get the index of the first nodal value associated with
	        // this FaceElement
	        unsigned first_index=
	         bnod_pt->index_of_first_value_assigned_by_face_element(this->Id);

	        // higher dimensional quantities
	        for(unsigned i=0;i<dim_el+1;i++)
	         {
	          x[i]+=this->nodal_position(j,i)*psi(j,0); // need to sort
	                                              // this out properly
	                                              // for generalised dofs
	          lambda[i]+=nod_pt->value
	           (first_index+i)*psi(j,0);
	         }
	        //In-element quantities
	        for(unsigned i=0;i<dim_el;i++)
	         {
	          //Loop over positional types
	          for (unsigned k=0;k<n_position_type;k++)
	           {
	            zeta[i]+=this->zeta_nodal(j,k,i)*psi(j,k);
	           }
	         }
	       }

	      // Get prescribed wall shape
	      Vector<double> r_prescribed(dim_el+1);
	      this->Boundary_shape_geom_object_pt->position(zeta,r_prescribed);


			}

				Vector<double> resultant_force(dim_el+1,0.0);
				this->compute_resultant_force(resultant_force);


			 	outfile << this->size() << " ";
			  outfile << resultant_force[0] << " ";
				outfile << resultant_force[1] << " ";
				outfile << resultant_force[2] << std::endl;


	   }


	 /// \short Compute the resultant force
 void compute_resultant_force(Vector<double> &resultant_force)
  {
   //Find out how many positional dofs there are
   unsigned n_position_type = this->nnodal_position_type();

#ifdef PARANOID
   if(n_position_type!=1)
    {
     throw OomphLibError(
      "ImposeDisplacementByLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   //Find out how many nodes there are
   unsigned n_node = this->nnode();

   // Dimension of element
   unsigned dim_el=this->dim();

   //Set up memory for the shape functions
   Shape psi(n_node);
   DShape dpsids(n_node,dim_el);

   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();


   // Initialise resultant force
   //Vector<double> resultant_force(dim_el+1,0.0);

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);

     //Only need to call the local derivatives
     this->dshape_local_at_knot(ipt,psi,dpsids);

     //Calculate the Eulerian coordinates and Lagrange multiplier
     Vector<double> x(dim_el+1,0.0);
     Vector<double> lambda(dim_el+1,0.0);
     Vector<double> zeta(dim_el,0.0);
     DenseMatrix<double> interpolated_a(dim_el,dim_el+1,0.0);

     // Loop over nodes
     for(unsigned j=0;j<n_node;j++)
      {
       Node* nod_pt=this->node_pt(j);

       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt =
        dynamic_cast<BoundaryNodeBase*>(this->node_pt(j));

       // Get the index of the first nodal value associated with
       // this FaceElement
       unsigned first_index=
        bnod_pt->index_of_first_value_assigned_by_face_element(this->Id);

       //Assemble higher-dimensional quantities
       for(unsigned i=0;i<dim_el+1;i++)
        {
         x[i]+=this->nodal_position(j,i)*psi(j);
         lambda[i]+=nod_pt->value(first_index+i)*psi(j);
         for(unsigned ii=0;ii<dim_el;ii++)
          {
           interpolated_a(ii,i) +=
            this->lagrangian_position(j,i)*dpsids(j,ii);
          }
        }
       if (!this->Sparsify)
        {
         for(unsigned k=0;k<n_position_type;k++)
          {
           //Assemble in-element quantities: boundary coordinate
           for(unsigned i=0;i<dim_el;i++)
            {
             zeta[i]+=this->zeta_nodal(j,k,i)*psi(j,k);
            }
          }
        }
      }

     if (this->Sparsify) zeta=this->Zeta_sub_geom_object[ipt];


     //Now find the local undeformed metric tensor from the tangent Vectors
     DenseMatrix<double> a(dim_el);
     for(unsigned i=0;i<dim_el;i++)
      {
       for(unsigned j=0;j<dim_el;j++)
        {
         //Initialise surface metric tensor to zero
         a(i,j) = 0.0;
         //Take the dot product
         for(unsigned k=0;k<dim_el+1;k++)
          {
           a(i,j) += interpolated_a(i,k)*interpolated_a(j,k);
          }
        }
      }


     //Find the determinant of the metric tensor
     double adet =0.0;
     switch(dim_el+1)
      {

      case 2:
       adet = a(0,0);
       break;

      case 3:
       adet = a(0,0)*a(1,1) - a(0,1)*a(1,0);
       break;

      default:
       throw
        OomphLibError(
         "Wrong dimension fill_in_generic_contribution_to_residuals_displ_lagr_multiplier",
         "ImposeDisplacementByLagrangeMultiplierElement::fill_in_generic_contribution_to_residuals_displ_lagr_multiplier()",
         OOMPH_EXCEPTION_LOCATION);
      }

     // Get prescribed wall shape
     Vector<double> r_prescribed(dim_el+1);
     if (!this->Sparsify)
      {
       this->Boundary_shape_geom_object_pt->position(zeta,r_prescribed);
      }
     else
      {
       this->Sub_geom_object_pt[ipt]->position(zeta,r_prescribed);
      }

     //Premultiply the weights and the square-root of the determinant of
     //the metric tensor
     double W = w*sqrt(adet);

     // Assemble error

     //Loop over directions
     for(unsigned i=0;i<dim_el+1;i++)
      {
       resultant_force[i] += -lambda[i]*W;
      }


    } //End of loop over the integration points


  }


};





  //===============================================================
  // A class for transversely isotropic constitutive laws -
  // inherits from isotropic constitutive class
  //===============================================================
  class TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw : public
  IsotropicStrainEnergyFunctionConstitutiveLaw
  {
  public:

    /// constructor takes a pointer to the strain energy function
    TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw(
							     StrainEnergyFunction* const &strain_energy_function_pt) :
      IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_function_pt) {}

    //===========================================================================
    /// Calculate the deviatoric part
    /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// Also return the contravariant deformed metric
    /// tensor and the determinant of the deformed metric tensor.
    /// Uses correct 3D invariants for 2D (plane strain) problems.
    /// This is the version for the pure incompressible formulation.
    //============================================================================
    void calculate_second_piola_kirchhoff_stress(
						 const DenseMatrix<double> &g,
						 const DenseMatrix<double> &G,
						 DenseMatrix<double> &sigma_dev,
						 DenseMatrix<double> &Gup,
						 double &detG,
						 Vector<double> &xi);

  };

//===============================================================
// member function which calculates the 2nd piola-kirchoff stress
//===============================================================
void TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw::
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g,
					const DenseMatrix<double> &G,DenseMatrix<double> &sigma_dev,
					DenseMatrix<double> &Gup, double &detG, Vector<double> &xi)
{
  // Error checking
#ifdef PARANOID
  error_checking_in_input(g,G,sigma_dev);
#endif

  // find the dimension of the problem
  unsigned dim = g.nrow();

#ifdef PARANOID
  if (dim==1 || dim==2)
    {
      std::string function_name =
	"TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw::";
      function_name += "calculate_second_piola_kirchhoff_stress()";

   throw OomphLibError(
		       "Check constitutive equations carefully when dim=1",
		       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

  // calculate the contravariant undeformed and deformed metric tensors
  // and get the determinants of the metric tensors
  DenseMatrix<double> gup(dim);
  //double  detg = calculate_contravariant(g,gup);
  (void)calculate_contravariant(g,gup);
  detG = calculate_contravariant(G,Gup);

  // Define the direction of the fibres in this transversely isotropic model
  // the fibres point in the direction of a vector M - alignment
  Vector<double> M(3,0.0);
  M[0] = 0.0;
  M[1] = 0.0;
  M[2] = 1.0;


  // calculate the strain invariants
  Vector<double> I(5,0.0);
  /// The third strain invariant is the volume change = detG/detg, set as 1
  I[2] = 1.0;
  // the first strain invariant is equal to g^{ij}G_{ij}
  // the second is equal to G^{rs}g_{rs}I3
  // the third is related to the direction of the fibres
  //  - all can be put within one set of loops
  for(unsigned i=0;i<dim;i++)
    {
      for(unsigned j=0;j<dim;j++)
	{
	  I[0] += gup(i,j)*G(i,j);
	  I[1] += Gup(i,j)*g(i,j);
	  I[3] += M[i]*M[j]*G(i,j);
	}
    }

  // the fifth strain invariant needs four loops
  // I5 = M_i * M_j * g^{rs}G_{is}G_{sj}
  for(unsigned i=0;i<dim;i++)
    {
      for(unsigned j=0;j<dim;j++)
	{
	  for(unsigned r=0;r<dim;r++)
	     {
	       for(unsigned s=0;s<dim;s++)
		 {
		   I[4] += M[i]*G(i,r)*gup(r,s)*G(s,j)*M[j];
		 }
	     }
	}
    }

  //Second strain invariant is multiplied by the third.
  //I[1] *= I[2];

  //Calculate the derivatives of the strain energy function wrt the
  //strain invariants
  Vector<double> dWdI(5,0.0);
   Strain_energy_function_pt->derivatives(I,dWdI);

	 //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
   //if the derivative wrt the second strain invariant is non-zero
   DenseMatrix<double> Bup(dim,dim,0.0);
   if(std::fabs(dWdI[1]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       Bup(i,j) = I[0]*gup(i,j);
	       for(unsigned r=0;r<dim;r++)
		 {
		   for(unsigned s=0;s<dim;s++)
		     {
		       Bup(i,j) -= gup(i,r)*gup(j,s)*G(r,s);
		     }
		 }
	     }
	 }
     }

   //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
   //if the derivative wrt the second strain invariant is non-zero
   DenseMatrix<double> Mup(dim,dim,0.0);
   if(std::fabs(dWdI[3]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       Mup(i,j) = M[i]*M[j];
	     }
	 }
     }

   //Only bother to compute the tensor C^{ij} (Green & Zerna notation)
   //if the derivative wrt the fifth strain invariant is non-zero
   DenseMatrix<double> Cup(dim,dim,0.0);
   if(std::fabs(dWdI[4]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       for(unsigned r=0;r<dim;r++)
		 {
		   for(unsigned s=0;s<dim;s++)
		     {
		       Cup(i,j) += M[i]*G(r,s)*gup(r,j)*M[s] + G(r,s)*gup(r,i)*M[s]*M[j];
		     }
		 }
	     }
	 }
     }

   // Now set the values of the functions phi and psi (Green and Zerna)
   // and chi and omega (derived due to anisotropy)
   double phi = 2.0*dWdI[0];
   double psi = 2.0*dWdI[1];
   double chi = 2.0*dWdI[3];
   double omega = 2.0*dWdI[4];

   double K;
   K = (I[0]*phi + 2.0*psi*I[1] + chi*I[3] + 2.0*omega*I[4])/3.0;

   // Now take the mechanical pressure P*G^{ij} away from sigma to
   // get the deviatoric part of the 2nd piola-kirchoff stress
   for(unsigned i=0;i<dim;i++)
     {
       for(unsigned j=0;j<dim;j++)
       {
	 sigma_dev(i,j) = phi*gup(i,j) + psi*Bup(i,j) + chi*Mup(i,j)
	   + omega*Cup(i,j) - K*Gup(i,j);
       }
     }


}


//===============================================================
// HGO (Holzapfel-Gasser-Ogden) strain energy function - this
// inherits from the standard strain energy function.
// W = (C1/2)*(I1-3) + k1/(2*k2)(exp(k2(I4-1)^2)-1)
//===============================================================
class HolzapfelGasserOgdenSEF : public StrainEnergyFunction
{

  public:

 /// Constructor takes the pointer to the value of the constants
 HolzapfelGasserOgdenSEF(double* c1_pt, double* k1_pt, double* k2_pt) : StrainEnergyFunction(),
   C1_pt(c1_pt), K1_pt(k1_pt), K2_pt(k2_pt)  {}


 /// Empty Virtual destructor
 virtual ~HolzapfelGasserOgdenSEF(){}

 /// Return the strain energy in terms of strain tensor
 double W(const DenseMatrix<double> &gamma)
  {return StrainEnergyFunction::W(gamma);}

 /// Return the strain energy in terms of the strain invariants
 double W(const Vector<double> &I)
  {return (*C1_pt)*(I[0]-3.0)/2.0 + ((*K1_pt)/(2*(*K2_pt)))*(exp((*K2_pt)*(pow(I[3]-1.0,2.0)))-1.0);}


 /// \short Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 void derivatives(Vector<double> &I, Vector<double> &dWdI)
  {
    dWdI[0] =(*C1_pt)/2.0;
   dWdI[1] = 0.0;
   dWdI[2] = 0.0;
   dWdI[3] = (*K1_pt)*(I[3]-1.0)*(exp((*K2_pt)*(pow(I[3]-1.0,2.0))));
   dWdI[4] = 0.0;
	 //std::cout << "lambda = " << sqrt(I[3]) << ", W = " << W(I) << std::endl;
  }

 /// \short Pure virtual function in which the user must declare if the
 /// constitutive equation requires an incompressible formulation
 /// in which the volume constraint is enforced explicitly.
 /// Used as a sanity check in PARANOID mode. True
  bool requires_incompressibility_constraint(){return true;}


  private:

 /// Pointer to constants
 double* C1_pt;
 double* K1_pt;
 double* K2_pt;

};



//====== Displaced Boundary ===============================
// Longitudinally displaced boundary in 3D
//===========================================================
class DisplacedBoundary : public GeomObject
{

public:

  /// constructor: specify amplitude of deflection from straight  horizontal line
  DisplacedBoundary(const double& ampl) : GeomObject(2,3)
  {
    Ampl=ampl;
  }

  /// Broken copy constructor
  DisplacedBoundary(const DisplacedBoundary& dummy)
  {
    BrokenCopy::broken_copy("DisplacedBoundary");
  }

  /// Broken assignment operator
  void operator=(const DisplacedBoundary&)
  {
    BrokenCopy::broken_assign("DisplacedBoundary");
  }

  /// Empty destructor
  ~DisplacedBoundary(){}


  // top surface of the cylinder shrinks
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {

	r[0] = zeta[0];
	r[1] = zeta[1];
	r[2] =  Global_Parameters::length + Ampl;

  }

  /// Short parameterised position on object: r(zeta). Evaluated at
  /// previous timestep.
  void position(const unsigned& t, const Vector<double>& zeta,
		Vector<double>& r) const
  {
    position(zeta,r);
  }

  /// Access to amplitude
  double& ampl() {return Ampl;}

  /// how many items of data does the shape of the object depend on?
  unsigned ngeom_data() const
  {
    return 0;
  }

private:

  /// amplitude of perturbation
  double Ampl;

};
}
