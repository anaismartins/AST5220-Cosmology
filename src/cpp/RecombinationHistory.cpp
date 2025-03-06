#include "RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo,
    double Yp, double z_reion, double delta_z_reion, double z_Hereion, double delta_z_Hereion) : cosmo(cosmo),
                                                                                                 Yp(Yp),
                                                                                                 z_reion(z_reion),
                                                                                                 delta_z_reion(delta_z_reion),
                                                                                                 z_Hereion(z_Hereion),
                                                                                                 delta_z_Hereion(delta_z_Hereion)

{
}

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
      // x_array.push_back(x_array[i]);
      Xe_arr.push_back(Xe_current);
      ne_arr.push_back(ne_current);
    }
    else
    {

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx)
      {
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      double Xe_ini = Xe_current;
      Vector Xe_ic{Xe_ini};

      Vector x_ode;

      // Check that i is greater than 0
      if (i > 0)
      {
        x_ode = {x_array[i - 1]};
        x_ode.push_back(x_array[i]);
      }
      else
      {
        printf("i is less than 0\n");
      }
      peebles_Xe_ode.solve(dXedx, x_ode, Xe_ic);

      // Verify that peebles_Xe_ode.solve correctly populates the internal data structure.
      printf("peebles_Xe_ode.get_data_by_component(0).size() = %d\n", peebles_Xe_ode.get_data_by_component(0).size());

      auto Xe_array = peebles_Xe_ode.get_data_by_component(0);

      // Ensure that Xe_array is not empty before accessing its elements.
      if (Xe_array.empty())
      {
        printf("Xe_array is empty\n");
      }

      double z = cosmo->get_z(x_array[i]);

      // Introduce reionization
      double f_He = Yp / (4. * (1 - Yp));
      double y_reion = pow(1. + z_reion, 1.5);
      double y = pow(1. + z, 1.5);
      double delta_y_reion = 3. / 2. * sqrt(1. + z_reion) * delta_z_reion;

      double Xe = Xe_array.back() + (1. + f_He) / 2. * (1 + tanh(y_reion - y) / delta_y_reion);

      // Introduce Helium reionization
      Xe = Xe + f_He / 2. * (1 + tanh((z_Hereion - z) / delta_z_Hereion));

      // Store the result
      Xe_arr.push_back(Xe);
    }
  }

  Vector log_Xe = log(Xe_arr);

  log_Xe_of_x_spline.create(x_array, log_Xe, "Xe");

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
  const double xhi0 = Constants.xhi0;
  const double xhi1 = Constants.xhi1;

  // Fetch cosmological parameters
  const double OmegaB0 = cosmo->get_OmegaB(0.0);
  const double H0 = cosmo->get_H0();
  const double TCMB0 = cosmo->get_TCMB(0.0);

  // Derived parameters
  const double rhoc0 = 3.0 * H0 * H0 / (8.0 * M_PI * G);

  double n_H = OmegaB0 * rhoc0 / (m_H * a * a * a);
  double nb = n_H / (1 - Yp);

  double fe0 = 1.0;

  double T_b = TCMB0 / a;
  double cHp = pow(m_e * T_b * Constants.k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-epsilon_0 / Constants.eV / (k_b * T_b)); // 1/m3
  double xHp = cHp / (fe0 * nb + cHp);

  double cHep = 2.0 * pow(m_e * T_b * Constants.k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-xhi0 / Constants.eV / (k_b * T_b));  // 1/m3
  double cHepp = 4.0 * pow(m_e * T_b * Constants.k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-xhi1 / Constants.eV / (k_b * T_b)); // 1/m3

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

  printf("x = %f, Xe = %f, ne = %f\n", x, Xe, ne);
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
  const double H_0 = cosmo->get_H0() / 1e5 * Constants.Mpc; // 1/s
  const double H = cosmo->H_of_x(x) / 1e5 * Constants.Mpc;  // 1/s
  const double TCMB0 = cosmo->get_TCMB(0.0);

  const double alpha = 1. / 137.0359992;

  double n_b = 3. * H_0 * H_0 / (8. * M_PI * G * m_H) * pow(a, -3) / 1e10 * Constants.Mpc * Constants.Mpc * pow(Constants.hbar, -3) * Constants.c * Constants.k_b * Constants.k_b; // 1/m^3
  double n_H = (1. - Yp) * n_b;                                                                                                                                                    // 1/m3
  double n_1s = (1. - X_e) * n_H;
  double lambda_alpha = H * pow(3. * epsilon_0, 3.) / (8. * M_PI * 8 * M_PI * n_1s) / Constants.eV * sqrt(pow(Constants.hbar / Constants.k_b, 5) / Constants.c); // 1/s

  double T_b = TCMB0 / a;
  double phi_2 = 0.448 * log(epsilon_0 / Constants.eV / (k_b * T_b));                                                                                                  // dimensionless
  double alpha_2 = 64. * M_PI / sqrt(27 * M_PI) * alpha * alpha / (m_e * m_e) * sqrt(epsilon_0 / (T_b * k_b)) * phi_2 * Constants.hbar * Constants.hbar / Constants.c; // m3/2
  double beta = alpha_2 * pow(m_e * T_b * Constants.k_b / (2. * M_PI * Constants.hbar * Constants.hbar), 1.5) * exp(-epsilon_0 / Constants.eV / (k_b * T_b));          // 1/s
  double beta2 = beta * exp(3. * epsilon_0 / (4. * k_b * T_b));
  double C_r = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2);

  dXedx[0] = C_r / H * (beta * (1 - X_e) - n_H * alpha_2 * X_e * X_e);

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
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;

    const double H = cosmo->H_of_x(x) / 1e5 * Constants.Mpc; // 1/s;

    double n_e = ne_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -c * n_e * sigma_T / H;

    return GSL_SUCCESS;
  };

  // Create a backwards x_array to integrate backwards in time
  Vector x_array_reversed(x_array.rbegin(), x_array.rend());

  double tau_ini = 0.0;
  Vector tau_ic{tau_ini};

  ODESolver ode;
  ode.solve(dtaudx, x_array_reversed, tau_ic);

  auto tau_array = ode.get_data_by_component(0);

  tau_of_x_spline.create(x_array_reversed, tau_array, "tau");

  // Compute visibility function
  Vector g_tilde_array(npts);
  for (int i = 0; i < npts; i++)
  {
    double x = x_array[i];
    double tau = tau_of_x(x);
    double g_tilde = exp(-tau) * dtaudx_of_x(x);
    g_tilde_array[i] = g_tilde;
  }

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
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const
{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const
{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const
{

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const
{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const
{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const
{
  double a = exp(x);

  const double G = Constants.G;
  const double m_H = Constants.m_H;

  const double H_0 = cosmo->get_H0() / 1e3 * Constants.Mpc; // 1/s

  double n_b = 3. * H_0 * H_0 / (8. * M_PI * G * m_H) * pow(a, -3) / 1e10 * Constants.Mpc * Constants.Mpc * pow(Constants.hbar, -3) * Constants.c * Constants.k_b * Constants.k_b; // 1/m^3
  double n_H = (1. - Yp) * n_b;                                                                                                                                                    // 1/m3

  return exp(log_Xe_of_x_spline(x)) * n_H;
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
