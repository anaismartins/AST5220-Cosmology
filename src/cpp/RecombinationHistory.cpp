#include "RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo,
    double Yp) : cosmo(cosmo),
                 Yp(Yp)
{
}

//====================================================
// Do all the solving we need to do - TODO: add reionization
//====================================================

void RecombinationHistory::solve()
{

  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons()
{
  Utils::StartTiming("Xe");

  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array;
  Vector Xe_arr;
  Vector ne_arr;

  // Calculate recombination history
  bool saha_regime = true;
  for (int i = 0; i < npts_rec_arrays; i++)
  {
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if (Xe_current < Xe_saha_limit)
      saha_regime = false;

    if (saha_regime)
    {
      x_array.push_back(x_array[i]);
      Xe_arr.push_back(Xe_current);
      ne_arr.push_back(ne_current);
    }
    else
    {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx)
      {
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result
      //=============================================================================
      //...
      //...
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x
  // functions are working
  //=============================================================================
  //...
  //...

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double, double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const
{
  const double a = exp(x);

  // Physical constants
  const double k_b = Constants.k_b;
  const double G = Constants.G;
  const double m_e = Constants.m_e;
  const double hbar = Constants.hbar;
  const double m_H = Constants.m_H;
  const double epsilon_0 = Constants.epsilon_0;
  const double H0_over_h = Constants.H0_over_h;

  const double chi0 = 24.5874 * Constants.eV;
  const double chi1 = 54.42279 * Constants.eV;

  // Fetch cosmological parameters
  const double OmegaB0 = cosmo->get_OmegaB(0.0);
  const double H0 = cosmo->get_H0();
  const double TCMB0 = cosmo->get_TCMB(0.0);

  // Derived parameters
  const double rhoc0 = 3.0 * H0 * H0 / (8.0 * M_PI * G);

  double nH = OmegaB0 * rhoc0 / (m_H * a * a * a);
  double nb = nH / (1 - Yp);

  // TODO: check units
  double fe0 = 1.0;

  double Tb = TCMB0 / a;
  double cHp = pow(m_e * Tb / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-epsilon_0 / (k_b * Tb));
  double xHp = cHp / (fe0 * nb + cHp);

  double cHep = 2.0 * pow(m_e * Tb / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-chi0 / (k_b * Tb));
  double cHepp = 4.0 * pow(m_e * Tb / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-chi1 / (k_b * Tb));

  double xHep = cHep * fe0 * nb / (fe0 * nb * fe0 * nb + fe0 * nb * cHep + cHep * cHepp);
  double xHepp = cHepp * xHep / (fe0 * nb);

  // double ne = 2.0 * nHepp + nHep + nHp;
  double fe = (2.0 * xHepp + xHep) * Yp / 4.0 + xHp * (1.0 - Yp);

  while ((abs(fe - fe0) > 1e-10))
  {
    fe0 = fe;
    xHp = cHp / (fe * nb + cHp);
    xHep = cHep * fe * nb / (fe * nb * fe * nb + fe * nb * cHep + cHep * cHepp);
    xHepp = cHepp * xHep / (fe * nb);
    fe = (2.0 * xHepp + xHep) * Yp / 4.0 + xHp * (1.0 - Yp);
  }

  double Xe = fe / (1 - Yp);
  double ne = nb * Xe;

  return std::pair<double, double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx)
{

  // Current value of a and X_e
  const double X_e = Xe[0];
  const double a = exp(x);

  // Physical constants in SI units
  const double k_b = Constants.k_b;
  const double G = Constants.G;
  const double c = Constants.c;
  const double m_e = Constants.m_e;
  const double hbar = Constants.hbar;
  const double m_H = Constants.m_H;
  const double sigma_T = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0 = Constants.epsilon_0;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...

  dXedx[0] = 0.0;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau()
{
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx)
  {
    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const
{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const
{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ne_of_x(double x) const
{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::get_Yp() const
{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const
{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp << "\n";
  std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const
{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  const double x_min = x_start;
  const double x_max = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&](const double x)
  {
    fp << x << " ";
    fp << Xe_of_x(x) << " ";
    fp << ne_of_x(x) << " ";
    fp << tau_of_x(x) << " ";
    fp << dtaudx_of_x(x) << " ";
    fp << ddtauddx_of_x(x) << " ";
    fp << g_tilde_of_x(x) << " ";
    fp << dgdx_tilde_of_x(x) << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
