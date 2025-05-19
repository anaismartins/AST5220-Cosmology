#include "PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo,
    RecombinationHistory *rec,
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : cosmo(cosmo),
                         rec(rec),
                         pert(pert),
                         A_s(A_s),
                         n_s(n_s),
                         kpivot_mpc(kpivot_mpc)
{
}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve()
{
  generate_bessel_function_splines();

  line_of_sight_integration();

  printf("line_of_sight_integration done\n");

  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");

  // printf("cell_TT_spline created\n");

  auto cell_TE = solve_for_cell(thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");

  // printf("cell_TE_spline created\n");

  auto cell_EE = solve_for_cell(thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines()
{
  Utils::StartTiming("besselspline");

  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());

  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for (size_t i = 0; i < ells.size(); i++)
  {
    const int ell = ells[i];

    double xmin_jell = 0.;
    double xmax_jell = 3000.;
    double deltax_jell = 2. * M_PI / 16.;
    Vector x_jell = Utils::linspace(xmin_jell, xmax_jell, int((xmax_jell - xmin_jell) / deltax_jell));
    Vector j_ell = Vector(x_jell.size());

    for (int xi = 0; xi < x_jell.size(); xi++)
    {
      // Compute the bessel function
      double j = Utils::j_ell(ell, x_jell[xi]);

      // std::cout << "j_ell[" << ell << "](" << x_jell[xi] << ") = " << j << std::endl;
      // Store the result in the vector
      j_ell[xi] = j;
    }

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(x_jell, j_ell, "j_ell_spline");

    // printf("j_ell_spline[%d] created\n", ell);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    const Vector &k_array,
    std::function<double(double, double)> &source_function)
{
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double nx = 300.;
  double xmin_los = -10.;
  double xmax_los = 0.;
  double deltax_los = (xmax_los - xmin_los) / nx;
  Vector x_los = Utils::linspace(xmin_los, xmax_los, int(nx));

  double eta0 = cosmo->eta_of_x(0.0);

  for (size_t ik = 0; ik < k_array.size(); ik++)
  {

    for (int elli = 0; elli < ells.size(); elli++)
    {
      const int ell = ells[elli];
      double integral = 0;

      for (int xi = 0; xi < x_los.size(); xi++)
      {
        // Compute the bessel function
        double eta = cosmo->eta_of_x(x_los[xi]);

        double jell = j_ell_splines[elli](k_array[ik] * (eta0 - eta));

        // Compute the source function
        double S = source_function(x_los[xi], k_array[ik]);

        // Do the integration
        double integral_x = S * jell * deltax_los;
        integral += integral_x;
      }

      // Store the result in result[ell][ik]
      result[elli][ik] = integral;
    }
  }

  printf("line_of_sight_integration_single done\n");

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration()
{
  const int nells = ells.size();

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double, double)> source_function_T = [&](double x, double k)
  {
    return pert->get_Source_T(x, k);
  };

  double eta0 = cosmo->eta_of_x(0.0);
  double deltak = 2. * M_PI / eta0 / 6.;
  double k_max = 3000. / eta0;
  int n_k = int((k_max - k_min) / deltak);

  Vector k_array = Utils::linspace(k_min, k_max, n_k);

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int elli = 0; elli < nells; elli++)
  {
    const int ell = ells[elli];
    thetaT_ell_of_k_spline[elli].create(k_array, thetaT_ell_of_k[elli], "Theta_ell(k)");
    // std::cout << "thetaT_ell_of_k_spline[" << ell << "] created\n";
    // printf("Theta_ell_spline[%d] created\n", ell);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if (Constants.polarization)
  {

    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    std::function<double(double, double)> source_function_E = [&](double x, double k)
    {
      return pert->get_Source_E(x, k);
    };

    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

    for (int elli = 0; elli < nells; elli++)
    {
      const int ell = ells[elli];
      thetaE_ell_of_k_spline[elli].create(k_array, thetaE_ell_of_k[elli], "ThetaE_ell(k)");
      // printf("ThetaE_ell_spline[%d] created\n", ell);
    }
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE)
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    std::vector<Spline> &f_ell_spline,
    std::vector<Spline> &g_ell_spline)
{
  const int nells = ells.size();
  Vector result = Vector(nells, 0.0);

  double factor = 4.0 * M_PI;

  double eta0 = cosmo->eta_of_x(0.0);
  double k_max = 3000. / eta0;

  double deltak_log = 2. * M_PI / (eta0 * k_max * 6.);
  int n_k_log = int((log(k_max) - log(k_min)) / deltak_log);

  Vector log_k_array = Utils::logspace(log10(k_min), log10(k_max), n_k_log);

  for (int elli = 0; elli < nells; elli++)
  {
    const int ell = ells[elli];

    for (int ik = 0; ik < log_k_array.size(); ik++)
    {
      double k = log_k_array[ik];
      double f_ell = f_ell_spline[elli](k);
      double g_ell = g_ell_spline[elli](k);

      // check for nans or other problematic values
      if (std::isnan(f_ell) || std::isnan(g_ell))
      {
        printf("f_ell[%d](%.2f) = %.2f\n", ell, k, f_ell);
        printf("g_ell[%d](%.2f) = %.2f\n", ell, k, g_ell);
        exit(1);
      }

      double primordial = primordial_power_spectrum(k);

      double integrand = primordial * f_ell * g_ell * deltak_log;

      result[elli] += integrand;
    }
    result[elli] *= factor;
  }

  printf("solve_for_cell done\n");
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const
{
  return A_s * pow(Constants.Mpc * k / kpivot_mpc, n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const
{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  double Phi = pert->get_Phi(x, k_mpc);
  double Omega_M = cosmo->get_OmegaB(x) + cosmo->get_OmegaCDM(x);
  double a = exp(x);
  double H0 = cosmo->H_of_x(0.0);

  double Delta_M = Constants.c * Constants.c * k_mpc * k_mpc * Phi / (3. / 2. * Omega_M * pow(a, -1.) * H0 * H0);
  double primordial = primordial_power_spectrum(k_mpc) / (k_mpc * k_mpc * k_mpc) * 2. * M_PI * M_PI;

  pofk = Delta_M * Delta_M * primordial;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const
{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const
{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const
{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const
{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size() - 1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax - 1);
  auto print_data = [&](const double ell)
  {
    double normfactor = (ell * (ell + 1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell + 1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB() * pow(4.0 / 11.0, 1.0 / 3.0), 2);
    double normfactorL = (ell * (ell + 1)) * (ell * (ell + 1)) / (2.0 * M_PI);
    fp << ell << " ";
    fp << thetaT_ell_of_k_spline[0](ell) * normfactor << " ";
    fp << cell_TT_spline(ell) * normfactor << " ";
    if (Constants.polarization)
    {
      fp << cell_EE_spline(ell) * normfactor << " ";
      fp << cell_TE_spline(ell) * normfactor << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

extern "C"
{
  PowerSpectrum *PowerSpectrum_new(
      BackgroundCosmology *cosmo,
      RecombinationHistory *rec,
      Perturbations *pert,
      double A_s,
      double n_s,
      double kpivot_mpc)
  {
    return new PowerSpectrum(cosmo, rec, pert, A_s, n_s, kpivot_mpc);
  }
  void PowerSpectrum_solve(PowerSpectrum *ps)
  {
    ps->solve();
  }
  double PowerSpectrum_get_thetaT_ell_of_k(
      PowerSpectrum *ps,
      const double ell,
      const double k)
  {
    // get index of ell
    int ell_index = -1;
    for (int i = 0; i < ps->ells.size(); i++)
    {
      if (ps->ells[i] == ell)
      {
        ell_index = i;
        break;
      }
    }
    return ps->thetaT_ell_of_k_spline[ell_index](k);
  }
  double PowerSpectrum_get_cell_TT(
      PowerSpectrum *ps,
      const double ell)
  {
    return ps->get_cell_TT(ell);
  }
  double PowerSpectrum_get_matter_power_spectrum(
      PowerSpectrum *ps,
      const double x,
      const double k_mpc)
  {
    return ps->get_matter_power_spectrum(x, k_mpc);
  }
}