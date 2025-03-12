#include "BackgroundCosmology.h"
#include "Utils.h"
#include <math.h>

//====================================================
// Constructors
//====================================================
BackgroundCosmology::BackgroundCosmology(
    double h,
    double OmegaB0,
    double OmegaCDM0,
    double OmegaK0,
    double Neff,
    double TCMB) : h(h),
                   OmegaB0(OmegaB0),
                   OmegaCDM0(OmegaCDM0),
                   OmegaK0(OmegaK0),
                   Neff(Neff),
                   TCMB(TCMB)
{
  H0 = Constants.H0_over_h * h;

  OmegaGamma0 = 2. * M_PI * M_PI / 30. * pow(Constants.k_b * TCMB, 4) / (pow(Constants.hbar, 3) * pow(Constants.c, 5)) * 8 * M_PI * Constants.G / (3 * H0 * H0);
  OmegaNu0 = Neff * 7. / 8. * pow(4. / 11., 4. / 3.) * OmegaGamma0;
  OmegaLambda0 = 1. - OmegaB0 - OmegaK0 - OmegaCDM0 - OmegaGamma0 - OmegaNu0;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve()
{

  Vector x_array;
  int npts = 1000;

  x_array = Utils::linspace(x_start, x_end, npts);

  ODESolver ode;

  Utils::StartTiming("Eta");

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx)
  {
    detadx[0] = Constants.c / Hp_of_x(x);

    return GSL_SUCCESS;
  };

  double etaini = 0.0;
  Vector eta_ic{etaini};

  ode.solve(detadx, x_array, eta_ic);

  auto eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta");

  Utils::EndTiming("Eta");

  Utils::StartTiming("CosmicTime");

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx)
  {
    dtdx[0] = 1. / H_of_x(x);

    return GSL_SUCCESS;
  };

  double tini = 1. / (2 * H_of_x(x_start));
  Vector t_ic{tini};

  ode.solve(dtdx, x_array, t_ic);

  auto t_array = ode.get_data_by_component(0);

  cosmic_time_spline.create(x_array, t_array, "cosmic time");

  Utils::EndTiming("CosmicTime");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const
{
  double a = exp(x);
  double H = H0 * sqrt((OmegaB0 + OmegaCDM0) * pow(a, -3) + (OmegaGamma0 + OmegaNu0) * pow(a, -4) + OmegaK0 * pow(a, -2) + OmegaLambda0);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const
{

  double a = exp(x);
  double H = H_of_x(x);

  return H * a;
}

double BackgroundCosmology::dHpdx_of_x(double x) const
{
  double a = exp(x);
  double OmegaM0 = OmegaB0 + OmegaCDM0;
  double OmegaR0 = OmegaGamma0 + OmegaNu0;

  double Hp = Hp_of_x(x);

  double dHpdx = Hp + 1. / 2. * (a * H0) * (a * H0) / Hp * (-3 * OmegaM0 * pow(a, -3) - 4 * OmegaR0 * pow(a, -4) - 2 * OmegaK0 * pow(a, -2));

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const
{
  double a = exp(x);
  double OmegaM0 = OmegaB0 + OmegaCDM0;
  double OmegaR0 = OmegaGamma0 + OmegaNu0;

  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);

  double ddHpddx = dHpdx * (1 - 1. / 4. * pow((a * H0) / Hp, 5) * pow((3 * OmegaM0 * pow(a, -3) + 4 * OmegaR0 * pow(a, -4) + 2 * OmegaK0 * pow(a, -2)), 2) * (9 * OmegaM0 * pow(a, -3) + 16 * OmegaR0 * pow(a, -4) + 4 * OmegaK0 * pow(a, -2)));

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const
{
  if (x == 0.0)
    return OmegaB0;

  double a = exp(x);
  double OmegaB = OmegaB0 / (pow(a, 3) * pow(H_of_x(x) / H0, 2));

  return OmegaB;
}

double BackgroundCosmology::get_OmegaGamma(double x) const
{
  if (x == 0.0)
    return OmegaGamma0;

  double a = exp(x);
  double OmegaGamma = OmegaGamma0 / (pow(a, 4) * pow(H_of_x(x) / H0, 2));
  return OmegaGamma;
}

double BackgroundCosmology::get_OmegaNu(double x) const
{

  if (x == 0.0)
    return OmegaNu0;

  double a = exp(x);
  double OmegaNu = OmegaNu0 / (pow(a, 4) * pow(H_of_x(x) / H0, 2));
  return OmegaNu;
}

double BackgroundCosmology::get_OmegaCDM(double x) const
{
  if (x == 0.0)
    return OmegaCDM0;

  double a = exp(x);
  double OmegaCDM = OmegaCDM0 / (pow(a, 3) * pow(H_of_x(x) / H0, 2));
  return OmegaCDM;
}

double BackgroundCosmology::get_OmegaLambda(double x) const
{
  if (x == 0.0)
    return OmegaLambda0;

  double a = exp(x);
  double OmegaLambda = OmegaLambda0 / (pow(H_of_x(x) / H0, 2));
  return OmegaLambda;
}

double BackgroundCosmology::get_OmegaK(double x) const
{
  if (x == 0.0)
    return OmegaK0;

  double a = exp(x);
  double OmegaK = OmegaK0 / (pow(a, 2) * pow(H_of_x(x) / H0, 2));
  return OmegaK;
}

double BackgroundCosmology::get_luminosity_distance_of_x(double x) const
{
  double a = exp(x);

  double chi = get_comoving_distance_of_x(x);
  double term = sqrt(abs(OmegaK0) * H0 * chi / Constants.c);

  double r;

  if (OmegaK0 == 0)
  {
    r = chi;
  }
  else if (OmegaK0 > 0)
  {
    r = chi * sinh(term) / term;
  }
  else
  {
    r = chi * sin(term) / term;
  }

  return r / a;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const
{
  double eta0 = eta_of_x(0.0);
  double etax = eta_of_x(x);

  return eta0 - etax;
}

double BackgroundCosmology::eta_of_x(double x) const
{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_cosmic_time(double x) const
{
  return cosmic_time_spline(x); // seconds
}

double BackgroundCosmology::get_z(double x) const
{
  return exp(x_end) / exp(x) - 1;
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
  std::cout << "OmegaB0:      " << OmegaB0 << "\n";
  std::cout << "OmegaCDM0:    " << OmegaCDM0 << "\n";
  std::cout << "OmegaLambda0: " << OmegaLambda0 << "\n";
  std::cout << "OmegaK0:      " << OmegaK0 << "\n";
  std::cout << "OmegaNu0:     " << OmegaNu0 << "\n";
  std::cout << "OmegaGamma0:  " << OmegaGamma0 << "\n";
  std::cout << "Neff:         " << Neff << "\n";
  std::cout << "h:            " << h << "\n";
  std::cout << "H0:           " << H0 << "\n";
  std::cout << "TCMB:         " << TCMB << "\n";
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
    fp << ddHpddx_of_x(x) << " ";
    fp << get_OmegaB(x) << " ";
    fp << get_OmegaCDM(x) << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaGamma(x) << " ";
    fp << get_OmegaNu(x) << " ";
    fp << get_OmegaK(x) << " ";
    fp << get_cosmic_time(x) << " ";
    fp << get_z(x) << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

extern "C"
{
  BackgroundCosmology *BackgroundCosmology_new(double h, double OmegaB0, double OmegaCDM0, double OmegaK0, double Neff, double TCMB)
  {
    return new BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB);
  }

  void BackgroundCosmology_solve(BackgroundCosmology *bc)
  {
    bc->solve();
  }

  double BackgroundCosmology_get_cosmic_time(BackgroundCosmology *bc, double x)
  {
    return bc->get_cosmic_time(x);
  }

  double BackgroundCosmology_get_OmegaGamma(BackgroundCosmology *bc, double x)
  {
    return bc->get_OmegaGamma(x);
  }

  double BackgroundCosmology_get_OmegaNu(BackgroundCosmology *bc, double x)
  {
    return bc->get_OmegaNu(x);
  }

  double BackgroundCosmology_get_H_of_x(BackgroundCosmology *bc, double x)
  {
    return bc->H_of_x(x);
  }

  double BackgroundCosmology_get_dHpdx_of_x(BackgroundCosmology *bc, double x)
  {
    return bc->dHpdx_of_x(x);
  }

  double BackgroundCosmology_eta_of_x(BackgroundCosmology *bc, double x)
  {
    return bc->eta_of_x(x);
  }
}