#include "RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo,
    double Yp, double z_reion, double delta_z_reion, double z_Hereion, double delta_z_Hereion, bool reion) : cosmo(cosmo),
                                                                                                             Yp(Yp),
                                                                                                             z_reion(z_reion),
                                                                                                             delta_z_reion(delta_z_reion),
                                                                                                             z_Hereion(z_Hereion),
                                                                                                             delta_z_Hereion(delta_z_Hereion),
                                                                                                             reion(reion)

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

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  Vector Xe_Saha(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  double x_change = 0.0;
  for (int i = 0; i < npts_rec_arrays; i++)
  {
    printf("i = %d\n", i);

    double Xe_current;
    double ne_current;

    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    if (std::isnan(Xe_ne_data.first) || std::isnan(Xe_ne_data.second))
    {
      printf("Invalid Xe or ne from Saha equation at index %d\n", i);
      exit(1);
    }

    // Electron fraction and number density
    Xe_current = Xe_ne_data.first;
    ne_current = Xe_ne_data.second;

    Xe_Saha[i] = Xe_current;

    // Are we still in the Saha regime?
    if (Xe_current < Xe_saha_limit)
    {
      saha_regime = false;
      if (x_change == 0.0)
      {
        x_change = x_array[i];
      }
    }

    if (saha_regime)
    {
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
    }
    else
    {
      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;

      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx)
      {
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      Vector x_ode, Xe_ic;
      if (i > 0)
      {
        double Xe_ini = Xe_arr[i - 1];
        Xe_ic = {Xe_ini};

        x_ode = {x_array[i - 1], x_array[i]};
      }
      else
      {
        double Xe_ini = Xe_current;
        Xe_ic = {Xe_ini};

        x_ode = {x_array[i], x_array[i] + 1e-6};
      }
      peebles_Xe_ode.solve(dXedx, x_ode, Xe_ic);

      auto Xe_array_Peebles = peebles_Xe_ode.get_data_by_component(0);

      if (std::isnan(Xe_array_Peebles[0]) || std::isnan(Xe_array_Peebles[1]))
      {
        printf("Invalid ODE solution at index %d\n", i);
        exit(1);
      }

      Xe_arr[i] = Xe_array_Peebles[1];
    }
  }

  printf("Reionization flag: %s\n", reion ? "true" : "false");

  if (reion)
  {
    for (int i = 0; i < npts_rec_arrays; i++)
    {
      if (x_array[i] > x_change)
      {
        double z = cosmo->get_z(x_array[i]);

        // Introduce reionization
        double f_He = Yp / (4. * (1 - Yp));
        double y_reion = pow(1. + z_reion, 1.5);
        double y = exp(-3. * x_array[i] / 2.);
        double delta_y_reion = 3. / 2. * sqrt(1. + z_reion) * delta_z_reion;

        // double f_e = Xe_array_Peebles[1] / (1. - Yp);

        double dXe_H = (1. + f_He) / 2. * (1 + tanh((y_reion - y) / delta_y_reion));
        double dXe_He = f_He / 2. * (1 + tanh((z_Hereion - z) / delta_z_Hereion));
        double Xe = Xe_arr[i] + dXe_H;

        // Introduce Helium reionization
        Xe = Xe + dXe_He;

        // std::cout << f_He << " " << z << " " << dXe_1 << " " << dXe_2 << " " << std::endl;

        // Store the result
        Xe_arr[i] = Xe;

        // check for nans
        if (std::isnan(Xe_arr[i]))
        {
          printf("Xe_arr[i] = %f\n", Xe_arr[i]);
          exit(0);
        }
      }
    }
  }

  Vector log_Xe = log(Xe_arr);

  // check for nans
  for (int i = 0; i < log_Xe.size(); i++)
  {
    if (std::isnan(log_Xe[i]))
    {
      printf("log_Xe[i] = %f\n", log_Xe[i]);
      exit(0);
    }
  }

  for (int i = 0; i < x_array.size(); i++)
  {
    if (std::isnan(x_array[i]) || std::isnan(log_Xe[i]))
    {
      printf("Invalid value in x_array or log_Xe at index %d\n", i);
      exit(1);
    }
  }

  log_Xe_of_x_spline.create(x_array, log_Xe, "Xe");
  Xe_of_x_spline_saha.create(x_array, Xe_Saha, "Xe_Saha");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double, double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const
{
  // increase print precision for the whole script
  std::cout.precision(15);

  // Current value of a and x
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

  if (std::isnan(Constants.k_b) || std::isnan(cosmo->get_H0()))
  {
    printf("Invalid constants or cosmological parameters\n");
    exit(1);
  }

  // Derived parameters
  const double rhoc0 = 3.0 * H0 * H0 / (8.0 * M_PI * G); // kg/m^3
  double n_H = OmegaB0 * rhoc0 / (m_H * a * a * a);      // 1/m^3
  double n_b = n_H / (1. - Yp);

  double fe = 1.0;
  double fe_old, n_e, cHp, cHep, cHepp, xHp, xHep, xHepp;

  double T_b = TCMB0 / a;

  double i = 0;
  do
  {
    fe_old = fe;
    n_e = fe * n_b;

    cHp = 1.0 / n_e * pow(m_e * T_b * k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-epsilon_0 / (k_b * T_b)); // 1/m3

    if (cHp < 1e-5)
    {
      double Xe = sqrt(cHp);
      double ne = n_b * Xe * (1 - Yp);
      return std::pair<double, double>(Xe, ne);
    }

    cHep = 2.0 / n_e * pow(m_e * T_b * k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-xhi0 / (k_b * T_b));  // 1/m3
    cHepp = 4.0 / n_e * pow(m_e * T_b * k_b / (2.0 * M_PI * hbar * hbar), 1.5) * exp(-xhi1 / (k_b * T_b)); // 1/m3

    xHp = cHp / (1.0 + cHp);
    xHep = cHep / (1.0 + cHep + cHep * cHepp);
    xHepp = cHepp * xHep;
    fe = (2.0 * xHepp + xHep) * Yp / 4.0 + xHp * (1.0 - Yp);

  } while ((abs(fe - fe_old) > 1e-10));

  double Xe = fe / (1. - Yp);
  double ne = n_b * Xe;

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
  const double H_0 = cosmo->get_H0(); // 1/s
  const double H = cosmo->H_of_x(x);  // 1/s
  const double TCMB0 = cosmo->get_TCMB(0.0);
  const double Omega_b0 = cosmo->get_OmegaB(0.0);

  const double alpha = 1. / 137.0359992;

  const double rho_c0 = 3. * H_0 * H_0 / (8. * M_PI * G);   // kg/m^3
  const double n_b = Omega_b0 * rho_c0 / (m_H * pow(a, 3)); // 1/m^3
  double n_H = (1. - Yp) * n_b;                             // 1/m3
  double n_1s = (1. - X_e) * n_H;
  double lambda_alpha = H * pow(3. * epsilon_0, 3.) / (8. * M_PI * 8 * M_PI * n_1s) / pow(hbar * c, 3); // 1/s

  double T_b = TCMB0 / a;
  double phi_2 = 0.448 * log(epsilon_0 / (k_b * T_b));                                                                                    // dimensionless
  double alpha_2 = 64. * M_PI / sqrt(27. * M_PI) * alpha * alpha / (m_e * m_e) * sqrt(epsilon_0 / (T_b * k_b)) * phi_2 * hbar * hbar / c; // m3/2

  double beta;
  if (epsilon_0 / (T_b * k_b) > 700)
  { // Prevent underflow
    beta = 0;
  }
  else
  {
    beta = alpha_2 * pow(m_e * T_b * k_b / (2. * M_PI * hbar * hbar), 1.5) * exp(-epsilon_0 / (k_b * T_b)); // 1/s
  }

  double beta2;
  if (3. * epsilon_0 / (4. * k_b * T_b) > 700)
  { // Prevent overflow
    beta2 = 0;
  }
  else
  {
    beta2 = beta * exp(3. * epsilon_0 / (4. * k_b * T_b));
  }
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

  // Make a reversed x-array
  auto x_array_reversed = x_array;
  std::reverse(x_array_reversed.begin(), x_array_reversed.end());

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx)
  {
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;

    const double H = cosmo->H_of_x(x); // 1/s;

    double n_e = ne_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -c * n_e * sigma_T / H;

    return GSL_SUCCESS;
  };

  // Create a backwards x_array to integrate backwards in time
  // Vector x_array_reversed(x_array.rbegin(), x_array.rend());

  double tau_ini = 0.0;
  Vector tau_ic{tau_ini};

  ODESolver ode;
  // ode.solve(dtaudx, x_array_reversed, tau_ic);
  ode.solve(dtaudx, x_array_reversed, tau_ic);

  auto tau_array = ode.get_data_by_component(0);

  // Reverse the tau array - now tau_array is in same order as x_array
  std::reverse(tau_array.begin(), tau_array.end());

  for (int i = 0; i < x_array.size(); i++)
  {
    if (std::isnan(x_array[i]) || std::isnan(tau_array[i]))
    {
      printf("Invalid value in x_array or tau_array at index %d\n", i);
      exit(1);
    }
  }

  // tau_of_x_spline.create(x_array_reversed, tau_array, "tau");
  tau_of_x_spline.create(x_array, tau_array, "tau");

  printf("Made tau spline\n");

  // Compute visibility function
  Vector g_tilde_array(npts);
  for (int i = 0; i < npts; i++)
  {
    double x = x_array[i];
    double g_tilde = -tau_of_x_spline.deriv_x(x) * exp(-tau_of_x_spline(x));
    // check for problems with g_tilde
    if (std::isnan(g_tilde) || std::isinf(g_tilde))
    {
      printf("dtaudx = %d, tau = %d\n", tau_of_x_spline.deriv_x(x), tau_of_x_spline(x));
      printf("Invalid value in g_tilde_array at index %d: %d\n", i, g_tilde);
      exit(1);
    }
    g_tilde_array[i] = g_tilde;
  }

  // check for problems with g_tilde
  for (int i = 0; i < g_tilde_array.size(); i++)
  {
    if (std::isnan(g_tilde_array[i]) || std::isinf(g_tilde_array[i]))
    {
      printf("Invalid value in g_tilde_array at index %d\n", i);
      exit(1);
    }
  }

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g");

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

double RecombinationHistory::Xe_of_x_Saha(double x) const
{
  return Xe_of_x_spline_saha(x);
}

double RecombinationHistory::ne_of_x(double x) const
{
  double a = exp(x);

  const double G = Constants.G;
  const double m_H = Constants.m_H;

  const double H_0 = cosmo->get_H0(); // 1/s
  const double Omega_b0 = cosmo->get_OmegaB(0.0);

  const double rho_c0 = 3. * H_0 * H_0 / (8. * M_PI * G);   // kg/m^3
  const double n_b = Omega_b0 * rho_c0 / (m_H * pow(a, 3)); // 1/m^3
  const double n_H = (1. - Yp) * n_b;                       // 1/m3

  return Xe_of_x(x) * n_H;
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
    fp << Xe_of_x_Saha(x) << " ";
    fp << ne_of_x(x) << " ";
    fp << tau_of_x(x) << " ";
    fp << dtaudx_of_x(x) << " ";
    fp << ddtauddx_of_x(x) << " ";
    fp << g_tilde_of_x(x) << " ";
    fp << dgdx_tilde_of_x(x) << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    // std::cout << "Finished writing data at x = " << x << "\n";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  std::cout << "Finished writing data to file " << filename << "\n";
}

extern "C"
{
  RecombinationHistory *RecombinationHistory_new(BackgroundCosmology *cosmo, double Yp, double z_reion, double delta_z_reion, double z_Hereion, double delta_z_Hereion, bool reion)
  {
    return new RecombinationHistory(cosmo, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, reion);
  }

  void RecombinationHistory_solve(RecombinationHistory *rec)
  {
    rec->solve();
  }

  double RecombinationHistory_tau_of_x(RecombinationHistory *rec, double x)
  {
    return rec->tau_of_x(x);
  }

  double RecombinationHistory_g_tilde_of_x(RecombinationHistory *rec, double x)
  {
    return rec->g_tilde_of_x(x);
  }

  double RecombinationHistory_Xe_of_x(RecombinationHistory *rec, double x)
  {
    return rec->Xe_of_x(x);
  }
}