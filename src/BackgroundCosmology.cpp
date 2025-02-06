#include "BackgroundCosmology.h"
#include "Utils.h"
#include <math.h>

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h,
    double OmegaB,
    double OmegaCDM,
    double OmegaK,
    double Neff,
    double TCMB) : h(h),
                   OmegaB(OmegaB),
                   OmegaCDM(OmegaCDM),
                   OmegaK(OmegaK),
                   Neff(Neff),
                   TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================
  double H0 = Constants.H0_over_h * h;

  double OmegaNu = Neff * 7. / 8. * pow(4. / 11., 4. / 3.) * OmegaR;
  double OmegaRtot = OmegaR + OmegaNu;

  double OmegaLambda = 1. - OmegaB - OmegaK - OmegaCDM - OmegaRtot;

  double getOmegaGamma
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve()
{
  Utils::StartTiming("Eta");

  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  Vector x_array;

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx)
  {
    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x); // TODO: Change this to a place where it computes it for all
  double H = H0 * sqrt((OmegaB + OmegaCDM) * pow(a, -3) + (OmegaR + OmegaNu) * pow(a, -4) + OmegaK * pow(a, -2) + OmegaLambda);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::dHpdx_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaB(double x) const
{
  if (x == 0.0)
    return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaR(double x) const
{

  double OmegaR0 = 2. * M_PI * M_PI / 30. * pow(Constants.k_b * TCMB, 4) / (pow(Constants.hbar, 3) * pow(Constants.c, 5)) * 8 * M_PI * Constants.G / (3 * H0 * H0);
  double a = exp(x);
  if (x == 0.0)
    return OmegaR0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  return OmegaR0 / (pow(a, 4) * pow(H_of_x(x) / H0, 2)); // TODO: this is wrong, it should be H(a), check if I should write that or other way to go about it
}

double BackgroundCosmology::get_OmegaNu(double x) const
{
  if (x == 0.0)
    return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const
{
  if (x == 0.0)
    return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaLambda(double x) const
{
  if (x == 0.0)
    return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaK(double x) const
{
  if (x == 0.0)
    return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_luminosity_distance_of_x(double x) const
{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const
{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const
{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const
{
  return H0;
}

double BackgroundCosmology::get_h() const
{
  return h;
}

double BackgroundCosmology::get_Neff() const
{
  return Neff;
}

double BackgroundCosmology::get_TCMB(double x) const
{
  if (x == 0.0)
    return TCMB;
  return TCMB * exp(-x);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const
{
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK << "\n";
  std::cout << "OmegaNu:     " << OmegaNu << "\n";
  std::cout << "OmegaR:      " << OmegaR << "\n";
  std::cout << "Neff:        " << Neff << "\n";
  std::cout << "h:           " << h << "\n";
  std::cout << "TCMB:        " << TCMB << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const
{
  const double x_min = -10.0;
  const double x_max = 0.0;
  const int n_pts = 100;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&](const double x)
  {
    fp << x << " ";
    fp << eta_of_x(x) << " ";
    fp << Hp_of_x(x) << " ";
    fp << dHpdx_of_x(x) << " ";
    fp << get_OmegaB(x) << " ";
    fp << get_OmegaCDM(x) << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x) << " ";
    fp << get_OmegaNu(x) << " ";
    fp << get_OmegaK(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
