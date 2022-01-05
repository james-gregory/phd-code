
namespace Global_Parameters
{

  // Holzapfel-Gasser-Ogden
  double C1 = 1.0;
  double K1 = 10.0;
  double K2 = 10.0;

  double A0 = 0.1; // radius

  // number of layers in the initial mesh
  unsigned nlayer = 5;

  // length of tendon
  double length = 1;

  // no of times the bulk elements will be refined uniformly
  int no_uniform = 0;
  // no of times the top elements will be refined
  int no_boundary = 0;

  double step_size = 0.005 * length;

  double max_strain = 0.1 * length;

  int counter;

}
