// simple script for modelling a non-linear, transversely isotropic
// hyperelastic cylinder (e.g. tendon).
// uses the HGO strain energy function

//Generic routines
#include "generic.h"
// Edited source code to pass xi through as a parameter
#include "MySolid.h"
// edited so that all five strain invariants are looped over
#include "MyConstitutive.h"
// The mesh
#include "meshes/tube_mesh.h"
// a separate header file containing only material parameters
#include "parameters.h"
// a header file which includes derived classes and functions which can be used anywhere
#include "MyDerivedClasses.h"
// include header file with the elastic-rupture constants in

using namespace oomph;


// Global Objects
//=================================================================
namespace Global_Objects
{
  /// pointer to strain energy function
  StrainEnergyFunction* Strain_energy_function_pt;

  /// Pointer to constitutive law
  ConstitutiveLaw* Constitutive_law_pt;

  /// GeomObject specifying the shape of the boundary. (flat)
  //class DisplacedBoundary;//(const double& ampl);
  DisplacedBoundary Boundary_geom_object(0.0);

} // end of namespace


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//===============================================================
/// Tube upgraded to become a solid mesh
//===============================================================
template<class ELEMENT>
class ElasticTubeMesh : public virtual TubeMesh<ELEMENT>,
			public virtual SolidMesh
{

 public:

  /// short constructor
  ElasticTubeMesh(GeomObject* volume_pt,
		  const Vector<double>& centreline_limits,
		  const Vector<double>& theta_positions,
		  const Vector<double>& radius_box,
		  const unsigned& nlayer,
		  TimeStepper* time_stepper_pt=
		  &Mesh::Default_TimeStepper) :
    TubeMesh<ELEMENT>(volume_pt,centreline_limits,theta_positions,radius_box,
		      nlayer,time_stepper_pt),
   SolidMesh()

  {
    // Assign the initial lagrangian coordinates
    set_lagrangian_nodal_coordinates();

    Vector<double> zeta(2);
    /// Loop over the top boundary
    {
      unsigned b=2;
      unsigned n_node = this->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  SolidNode* nod_pt=static_cast<SolidNode*>(this->boundary_node_pt(b,n));
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
      this->Boundary_coordinate_exists[b]=true;
    }


  }


  /// Empty Destructor
  virtual ~ElasticTubeMesh() { }

};

//===============================================================
/// Tube upgraded to become a refineable solid mesh
//===============================================================
template<class ELEMENT>
class RefineableElasticTubeMesh : public virtual  RefineableTubeMesh<ELEMENT>,
				  public virtual  SolidMesh
{

 public:

  /// short constructor
  RefineableElasticTubeMesh(GeomObject* volume_pt,
			    const Vector<double>& centreline_limits,
			    const Vector<double>& theta_positions,
			    const Vector<double>& radius_box,
			    const unsigned& nlayer,
			    TimeStepper* time_stepper_pt=
			    &Mesh::Default_TimeStepper) :
    RefineableTubeMesh<ELEMENT>(volume_pt,centreline_limits,theta_positions,radius_box,
				nlayer,time_stepper_pt),
   SolidMesh()
  {

    this->setup_octree_forest();

    /// Assign the initial lagrangian coordinates
    set_lagrangian_nodal_coordinates();

    Vector<double> zeta(2);
    /// Loop over the top boundary
    {
      unsigned b=2;
      unsigned n_node = this->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  Node* nod_pt=this->boundary_node_pt(b,n);
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
      this->Boundary_coordinate_exists[b]=true;
    }

  }

  /// Empty Destructor
  virtual ~RefineableElasticTubeMesh() { }

};


//=start_of_MyCylinder===================================
//A geometric object that represents the geometry of the domain
//================================================================
class MyEllipticalCylinder : public GeomObject
{
public:

 /// Constructor that takes the radius of the tube as its argument
 MyEllipticalCylinder(const double &radius1, const double &radius2) :
  GeomObject(3,3), Radius1(radius1), Radius2(radius2) { }

/// Destructor
virtual~MyEllipticalCylinder(){}

///Lagrangian coordinate xi
void position (const Vector<double>& xi, Vector<double>& r) const
{
 r[0] = xi[2]*Radius1*cos(xi[1]);
 r[1] = xi[2]*Radius2*sin(xi[1]);
 r[2] = xi[0];}


/// Return the position of the tube as a function of time
/// (doesn't move as a function of time)
void position(const unsigned& t,
              const Vector<double>& xi, Vector<double>& r) const
  {
   position(xi,r);
  }

private:

 ///Storage for the radius of the tube
 double Radius1;
 double Radius2;

};






//==================================================================
// Deformation of elastic tube
//==================================================================
template<class ELEMENT>
class CylinderExtensionProblem : public Problem
{

public:

  /// Constructor:
  CylinderExtensionProblem();

  // Run simulation
  void run(const std::string &dirname);

  /// doc the solution
  void doc_solution(DocInfo& doc_info);

  /// update function (empty)
  void actions_after_newton_solve() {}

 #ifdef NOREFINE

  /// Access function for the solid mesh
  ElasticTubeMesh<ELEMENT>*& solid_mesh_pt()
  {return Solid_mesh_pt;}

 #else

  /// Access function for the refineable solid mesh
  RefineableElasticTubeMesh<ELEMENT>*& solid_mesh_pt()
  {return Solid_mesh_pt;}

 #endif

 #ifndef NOREFINE

  /// Actions before adapt
  void actions_before_adapt()
  {
    /// kill the elements and wipe the surface mesh
    delete_lagrange_multiplier_elements();
    /// Rebuild the problems global mesh from its various sub meshes
    rebuild_global_mesh();

  } // end of actions_before_adapt

  /// Actions after adapt
  void actions_after_adapt()
  {
	  {
	    unsigned b=0;
	    unsigned n_node = Solid_mesh_pt->nboundary_node(b);

				for(unsigned n=0;n<n_node;n++){
					// Pin nodes in z-direction and to prevent rotations
					// pin only in z direction and to prevent rotations and translations:
					//std::cout << "only pinning boundary nodes in z-direction" << std::endl;
					const double tolerance = 1e-10;
					// pin in the z-direction
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
					// if the node lies on the x axis, pin in y direction
					if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) < tolerance)
					&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) > (-tolerance)))
					{
						Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
					}
					// if the node lies on the y axis, pin in the x direction
					if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) <  tolerance)
					&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) > (-tolerance)))
					{
						Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
					}
				}

			}

	  // end of setup boundary conditions


    /// Loop over the top boundary
    {
      Vector<double> zeta(2);
      unsigned b=2;
      unsigned n_node = solid_mesh_pt()->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  SolidNode* nod_pt=static_cast<SolidNode*>(solid_mesh_pt()->boundary_node_pt(b,n));
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
    }

    /// Create the elements that impose the displaced boundary constraint
    /// and attach them to the bulk elements that are adjacent to boundary 2
    create_lagrange_multiplier_elements();

    /// Rebuild the problem's global mesh from its various sub-meshes
    rebuild_global_mesh();

    /// Pin the redundant solid pressures // PVDEquationsBase
    PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(Solid_mesh_pt->element_pt());


  } // end of  after adapt

 #endif


private:

  /// Create elements that enforce prescribed boundary motion
  void create_lagrange_multiplier_elements();

  /// delete elements that enforce prescribed boundary motion
  void delete_lagrange_multiplier_elements();

  /// updated before solve: empty
  void actions_before_newton_solve(){}


 #ifdef NOREFINE

  /// pointer to solid mesh
  ElasticTubeMesh<ELEMENT>* Solid_mesh_pt;

 #else

  /// pointer to solid mesh
  RefineableElasticTubeMesh<ELEMENT>* Solid_mesh_pt;

 #endif

  /// pointers to meshes of lagrange multiplier elements
  SolidMesh* Lagrange_multiplier_mesh_pt;

  /// pointer to GeomObject that specifies the domain volume
  GeomObject *Volume_pt;

  /// Docinfo object for output
  DocInfo Doc_info;

};

//==============================================================
// Creates Lagrange multiplier elements
//==============================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
  /// Lagrange multiplier elements are located on boundary 2
  unsigned b=2;

  // How many bulk elements are adjacent to boundary b?
  unsigned n_element = solid_mesh_pt()->nboundary_element(b);

  // loop over the bulk elements adjacent to boundary b
  for(unsigned e=0;e<n_element;e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary b
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->boundary_element_pt(b,e));

      // find the index of the face of element e along boundary b
      int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);


      /// create new element and add to mesh
      Lagrange_multiplier_mesh_pt->add_element_pt( new MyLagrangeMultiplierElement<ELEMENT>(bulk_elem_pt,face_index));

    }

  // Loop over the elements in the lagrange multiplier element mesh

  n_element = Lagrange_multiplier_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {

      ///cast to lagrange multiplier element
      MyLagrangeMultiplierElement<ELEMENT> *el_pt =
	dynamic_cast<MyLagrangeMultiplierElement<ELEMENT>*>
	(Lagrange_multiplier_mesh_pt->element_pt(i));

      /// Set the GeomObject that defines the boundary shape and
      /// specify which bulk boundary we are attached to ( need to extract
      /// the boundary coordinate from the bulk nodes)
      el_pt->set_boundary_shape_geom_object_pt(&Global_Objects::Boundary_geom_object,b);

		}




}// end of create_lagrane_multiplier_elements

///===============================================================
/// Deletes lagrange multiplier elements
///===============================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::delete_lagrange_multiplier_elements()
{
  /// How many surface elements are in the surface mesh?
  unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();

  /// Loop over the surface elements
  for(unsigned e=0;e<n_element;e++)
    {
      /// kill the surface element
      delete Lagrange_multiplier_mesh_pt->element_pt(e);
    }

  /// wipe the mesh
  Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();

} // end of delete_lagrange_multiplier_elements

///===============================================================
/// Constructor
///===============================================================
template<class ELEMENT>
CylinderExtensionProblem<ELEMENT>::CylinderExtensionProblem()
{
  /// Create GeomObject that specifies the domain geometry
  /// the radius of the reference cylinder

  Volume_pt = new MyEllipticalCylinder(Global_Parameters::A0, Global_Parameters::A0);

  /// Define pi
  const double pi = 4.0*atan(1.0);


  /// Set the centreline coordinates spanning the mesh
  Vector<double> centreline_limits(2);
  centreline_limits[0] = 0.0;
  centreline_limits[1] = Global_Parameters::length;

  /// Set the positions of the angles that divide the outer ring
  Vector<double> theta_positions(4);
  theta_positions[0] = -0.75*pi;
  theta_positions[1] = -0.25*pi;
  theta_positions[2] = 0.25*pi;
  theta_positions[3] = 0.75*pi;


  /// Define the radial fraction of the central box
  Vector<double> radial_frac(4,0.5);

#ifdef NOREFINE


  /// now create the mesh
  Solid_mesh_pt = new ElasticTubeMesh<ELEMENT>
    (Volume_pt,centreline_limits,theta_positions,radial_frac,Global_Parameters::nlayer);


#else
  /// now create the mesh
  Solid_mesh_pt = new RefineableElasticTubeMesh<ELEMENT>
    (Volume_pt,centreline_limits,theta_positions,radial_frac,Global_Parameters::nlayer);

  // /// Set error estimator
  Solid_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  Solid_mesh_pt->max_permitted_error() = 0.002;
  Solid_mesh_pt->min_permitted_error() = 0.0001;


#endif

  // set these positions as the unstrained node positions -- important when
  // applying transformations to the mesh
  Solid_mesh_pt->set_lagrangian_nodal_coordinates();
  /// loop over the elements in the mesh to set parameters/functions pointers
  unsigned n_element=Solid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      /// cast to a solid element
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));

      /// set the constitutive law
      el_pt->constitutive_law_pt() =
	Global_Objects::Constitutive_law_pt;

      /// Material is incompressible
      PVDEquationsWithPressure<3>* cast_el_pt =
	dynamic_cast<PVDEquationsWithPressure<3>*>(Solid_mesh_pt->element_pt(i));
      if (cast_el_pt!=0)
	{
	  cast_el_pt->set_incompressible();
	}

    }


  // REFINE ELEMENTS
  // refine all or just the top elements a specified number of times
  for(int i=0; i<Global_Parameters::no_uniform; i++)
    {
      Solid_mesh_pt->refine_uniformly();
    }

  std::cout << "number of elements: "<< Solid_mesh_pt->nelement() << std::endl;
	std::cout << "number of nodes: "<< Solid_mesh_pt->nnode() << std::endl;

  /// Construct the mesh of elements that enforce prescribed boundary motion
  Lagrange_multiplier_mesh_pt = new SolidMesh;
  create_lagrange_multiplier_elements();


  /// Solid mesh is the first sub-mesh
  add_sub_mesh(Solid_mesh_pt);

  /// add lagrange multiplier sub-mesh
  add_sub_mesh(Lagrange_multiplier_mesh_pt);

  /// build combined global mesh
  build_global_mesh();

  Solid_mesh_pt->disable_adaptation();
  /// setup boundary conditions
  {
		unsigned b=0;
		unsigned n_node = Solid_mesh_pt->nboundary_node(b);


			for(unsigned n=0;n<n_node;n++){
				// Pin nodes in z-direction and to prevent rotations
				// pin only in z direction and to prevent rotations and translations:
				//std::cout << "only pinning boundary nodes in z-direction" << std::endl;
				const double tolerance = 1e-10;
				// pin in the z-direction
				Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
				// if the node lies on the x axis, pin in y direction
				if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) < tolerance)
				&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) > (-tolerance)))
				{
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
				}
				// if the node lies on the y axis, pin in the x direction
				if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) <  tolerance)
				&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) > (-tolerance)))
				{
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
				}
			}

	} // end of setup boundary conditions

  /// pin the redundant solid pressures
  PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(Solid_mesh_pt->element_pt());

  /// attach the boundary conditions to the mesh
  cout << "Number of equations: " << assign_eqn_numbers() << std::endl;


}

///==========================================================
/// Doc the solution
///==========================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  ofstream some_file;
  char filename[150];

  /// number of plot points
  unsigned npts = 5;

  ///  the shape  of the deformed boday
  sprintf(filename, "%s/soln%s.dat", doc_info.directory().c_str(),
	  doc_info.label().c_str());
  some_file.open(filename);
  Solid_mesh_pt->output(some_file,npts);
  some_file.close();

  // output information about the surface elements on the top of the tendon
  sprintf(filename, "%s/resultant_force%s.dat",doc_info.directory().c_str(),
	  doc_info.label().c_str());
  some_file.open(filename);
  Lagrange_multiplier_mesh_pt->output(some_file,npts);
  some_file.close();



}


//======================================================
/// Driver for simple elastic problem
//======================================================
int main(int argc, char **argv)
{

	// store command line arguments
	CommandLineArgs::setup(argc,argv);
	// C1
	CommandLineArgs::specify_command_line_flag("--c1",&Global_Parameters::C1);
	// K1
	CommandLineArgs::specify_command_line_flag("--k1",&Global_Parameters::K1);
	// K2
	CommandLineArgs::specify_command_line_flag("--k2",&Global_Parameters::K2);
	// counter for the naming of files
	CommandLineArgs::specify_command_line_flag("--counter",&Global_Parameters::counter);
	// Parse command line
	CommandLineArgs::parse_and_assign();
	// Doc what has actually been specified on the command line
	CommandLineArgs::doc_specified_flags();

	// select the transversely isotropic HGO model
	Global_Objects::Strain_energy_function_pt =
	new HolzapfelGasserOgdenSEF(&Global_Parameters::C1,
	&Global_Parameters::K1,
	&Global_Parameters::K2);

	Global_Objects::Constitutive_law_pt =
	new TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw(
	Global_Objects::Strain_energy_function_pt);

	// create a file name and label using the ellipticity and taper
	string filename = "RESLT/data";
	std::cout << filename << std::endl;

	#ifdef NOREFINE
	///set up the problem with pressure/displacement formulation
	CylinderExtensionProblem<MyQPVDElementWithContinuousPressure<3> > problem;
	//problem.run(filename);

	#else
	/// set up the problem with pressure/displacement formulation
	CylinderExtensionProblem<MyRefineableQPVDElementWithContinuousPressure<3> > problem;
	//problem.run(filename);

	#endif

	// output
	DocInfo doc_info;
	// set output director
	doc_info.set_directory(filename);


	/// Initial parameter value
	Global_Objects::Boundary_geom_object.ampl()=0.0;


	/// doc initial configuration
	problem.doc_solution(doc_info);

	while(Global_Objects::Boundary_geom_object.ampl() < Global_Parameters::max_strain)
	{
		// increment the boundary displacement
		Global_Objects::Boundary_geom_object.ampl() += Global_Parameters::step_size;
		// solve
		problem.newton_solve(0);
		// document
		problem.doc_solution(doc_info);

	}




}
