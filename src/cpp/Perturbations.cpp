#include "Perturbations.h"
#include <algorithm> // For std::min

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo,
    RecombinationHistory *rec) : cosmo(cosmo),
                                 rec(rec)
{
}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve()
{

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations()
{
  Utils::StartTiming("integrateperturbation");

  Vector k_array(n_k);
  k_array = Utils::logspace(log10(k_min), log10(k_max), n_k);
  Vector x_array(n_x);
  x_array = Utils::linspace(x_start, x_end, n_x);

  Vector delta_cdm(n_x * n_k);
  Vector delta_b(n_x * n_k);
  Vector v_cdm(n_x * n_k);
  Vector v_b(n_x * n_k);
  Vector Phi(n_x * n_k);
  Vector Psi(n_x * n_k);
  Vector2D Theta = Vector2D(Constants.n_ell_theta, Vector(n_x * n_k));
  Vector2D Theta_p = Vector2D(Constants.n_ell_thetap, Vector(n_x * n_k));
  Vector Pi(n_x * n_k);
  Vector2D Nu = Vector2D(Constants.n_ell_neutrinos, Vector(n_x * n_k));

  // Loop over all wavenumbers
  for (int ik = 0; ik < n_k; ik++)
  {
    // Progress bar...
    if ((10 * ik) / n_k != (10 * ik + 10) / n_k)
    {
      std::cout << (100 * ik + 100) / n_k << "% " << std::flush;
      if (ik == n_k - 1)
        std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    // Integrate from x_start -> x_end_tight
    double i_tc;
    for (int ix = 0; ix < n_x; ix++)
    {
      if (x_array[ix] > x_end_tight)
      {
        i_tc = ix;
        break;
      }
    }

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx)
    {
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    ODESolver ode_tc;

    Vector x_array_tc(x_array.begin(), x_array.begin() + i_tc);
    Vector y_ic{y_tight_coupling_ini};

    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_ic);
    auto y_tight_coupling = ode_tc.get_data();

    for (int ix = 0; ix < i_tc; ix++)
    {
      double a = exp(x_array[ix]);

      delta_cdm[ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_deltacdm_tc];
      delta_b[ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_deltab_tc];
      v_cdm[ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_vcdm_tc];
      v_b[ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_vb_tc];
      Phi[ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_Phi_tc];
      for (int l = 0; l < Constants.n_ell_neutrinos_tc; l++)
      {
        Nu[l][ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_start_nu_tc + l];
      }
      for (int l = 0; l < Constants.n_ell_theta_tc; l++)
      {
        Theta[l][ix + n_x * ik] = y_tight_coupling[ix][Constants.ind_start_theta_tc + l];
      }
      Psi[ix + n_x * ik] = -Phi[ix + n_x * ik] - 12. * H0 * H0 / (c * c * k * k * a * a) * (Omega_nu0 * Nu[2][ix + n_x * ik]);

      Pi[ix + n_x * ik] = 0;
    }

    ODESolver ode_full;

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    Vector y_tc = ode_tc.get_data_by_xindex(i_tc - 1);

    auto y_full_ini = set_ic_after_tight_coupling(y_tc, x_end_tight, k);

    // check for nans in the initial conditions
    for (int i = 0; i < y_full_ini.size(); i++)
    {
      if (std::isnan(y_full_ini[i]))
      {
        printf("y_full_ini[%d] is nan\n", i);
        throw std::runtime_error("y_full_ini is nan");
      }
    }

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx)
    {
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_array_full(x_array.begin() + i_tc, x_array.end());
    Vector y_ic_full{y_full_ini};

    // check for nans
    for (int i = 0; i < y_ic_full.size(); i++)
    {
      if (std::isnan(y_ic_full[i]))
      {
        printf("y_ic_full[%d] is nan\n", i);
        throw std::runtime_error("y_ic_full is nan");
      }
    }

    ode_full.solve(dydx_full, x_array_full, y_ic_full);
    auto y_full = ode_full.get_data();

    // check y_full for nans
    for (int ix = 0; ix < y_full.size(); ix++)
    {
      if (std::isnan(y_full[ix][0]))
      {
        printf("y_full[%d][0] is nan\n", ix);
        throw std::runtime_error("y_full is nan");
      }
    }

    for (int ix = i_tc; ix < n_x; ix++)
    {
      double a = exp(x_array[ix]);

      delta_cdm[ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_deltacdm];
      delta_b[ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_deltab];
      v_cdm[ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_vcdm];
      v_b[ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_vb];
      Phi[ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_Phi];
      Psi[ix + n_x * ik] = -Phi[ix + n_x * ik] - 12. * H0 * H0 / (c * c * k * k * a * a) * (Omega_gamma0 * Theta[2][ix + n_x * ik] + Omega_nu0 * Nu[2][ix + n_x * ik]);
      for (int l = 0; l < Constants.n_ell_theta; l++)
      {
        Theta[l][ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_start_theta + l];
      }
      for (int l = 0; l < Constants.n_ell_thetap; l++)
      {
        Theta_p[l][ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_start_thetap + l];
      }
      Pi[ix + n_x * ik] = Theta[2][ix + n_x * ik] + Theta_p[0][ix + n_x * ik] + Theta_p[2][ix + n_x * ik];
      for (int l = 0; l < Constants.n_ell_neutrinos; l++)
      {
        Nu[l][ix + n_x * ik] = y_full[ix - i_tc][Constants.ind_start_nu + l];
      }
    }
  }

  Utils::EndTiming("integrateperturbation");

  delta_cdm_spline.create(x_array, k_array, delta_cdm, "delta_cdm_spline");
  delta_b_spline.create(x_array, k_array, delta_b, "delta_b_spline");
  v_cdm_spline.create(x_array, k_array, v_cdm, "v_cdm_spline");
  v_b_spline.create(x_array, k_array, v_b, "v_b_spline");
  Phi_spline.create(x_array, k_array, Phi, "Phi_spline");
  Psi_spline.create(x_array, k_array, Psi, "Psi_spline");
  Pi_spline.create(x_array, k_array, Pi, "Pi_spline");

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
  Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);

  // check for nans in theta
  for (int l = 0; l < Constants.n_ell_theta; l++)
  {
    for (int ix = 0; ix < n_x * n_k; ix++)
    {
      if (std::isnan(Theta[l][ix]))
      {
        printf("Theta[%d][%d] is nan\n", l, ix);
      }
    }
  }

  for (int l = 0; l < Constants.n_ell_theta; l++)
  {
    Theta_spline[l].create(x_array, k_array, Theta[l], "Theta_spline");
  }
  for (int l = 0; l < Constants.n_ell_thetap; l++)
  {
    Theta_p_spline[l].create(x_array, k_array, Theta_p[l], "Theta_p_spline");
  }
  for (int l = 0; l < Constants.n_ell_neutrinos; l++)
  {
    Nu_spline[l].create(x_array, k_array, Nu[l], "Nu_spline");
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const
{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc = Constants.n_ell_tot_tc;
  const bool polarization = Constants.polarization;
  const bool neutrinos = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm = y_tc[Constants.ind_deltacdm_tc];
  double &delta_b = y_tc[Constants.ind_deltab_tc];
  double &v_cdm = y_tc[Constants.ind_vcdm_tc];
  double &v_b = y_tc[Constants.ind_vb_tc];
  double &Phi = y_tc[Constants.ind_Phi_tc];
  double *Theta = &y_tc[Constants.ind_start_theta_tc];
  double *Nu = &y_tc[Constants.ind_start_nu_tc];

  double Hp = cosmo->Hp_of_x(x);

  double a = exp(x);
  double f_nu = Omega_nu0 / (Omega_gamma0 + Omega_nu0);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  const double Psi = -1. / (3. / 2. + 2. * f_nu / 5.);
  Phi = -(1. + 2. * f_nu / 5.) * Psi;
  delta_cdm = -3. / 2. * Psi;
  delta_b = delta_cdm;
  v_cdm = -c * k / (2. * Hp) * Psi;
  v_b = v_cdm;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -1. / 2. * Psi;
  Theta[1] = c * k / (6. * Hp) * Psi;

  // SET: Neutrino perturbations (N_ell)
  if (neutrinos)
  {
    Nu[0] = -1. / 2. * Psi;
    Nu[1] = c * k / (6. * Hp) * Psi;
    Nu[2] = -c * c * k * k * a * a * (Phi + Psi) / (12. * H0 * H0 * Omega_nu0);
    for (int l = 3; l < n_ell_neutrinos_tc; l++)
    {
      Nu[l] = c * k / ((2 * l + 1) * Hp) * Nu[l - 1];
    }
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc,
    const double x,
    const double k) const
{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta = Constants.n_ell_theta;
  const int n_ell_thetap = Constants.n_ell_thetap;
  const int n_ell_neutrinos = Constants.n_ell_neutrinos;
  const bool polarization = Constants.polarization;
  const bool neutrinos = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc = y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc = y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc = y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc = y_tc[Constants.ind_vb_tc];
  const double &Phi_tc = y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm = y[Constants.ind_deltacdm_tc];
  double &delta_b = y[Constants.ind_deltab_tc];
  double &v_cdm = y[Constants.ind_vcdm_tc];
  double &v_b = y[Constants.ind_vb_tc];
  double &Phi = y[Constants.ind_Phi_tc];
  double *Theta = &y[Constants.ind_start_theta_tc];
  double *Theta_p = &y[Constants.ind_start_thetap_tc];
  double *Nu = &y[Constants.ind_start_nu_tc];

  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  double f_nu = Omega_nu0 / (Omega_gamma0 + Omega_nu0);

  const double Psi = -1. / (3. / 2. + 2. * f_nu / 5.);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;

  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -8. * c * k / (15. * Hp * dtaudx) * Theta[1];
  for (int l = 3; l < n_ell_theta; l++)
  {
    Theta[l] = -l / (2. * l + 1.) * c * k / (Hp * dtaudx) * Theta[l - 1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if (polarization)
  {
    Theta_p[0] = 5. / 4. * Theta[2];
    Theta_p[1] = -c * k / (4. * Hp * dtaudx) * Theta[2];
    Theta_p[2] = 1. / 4. * Theta[2];
    for (int l = 3; l < n_ell_thetap; l++)
    {
      Theta_p[l] = -l / (2 * l + 1) * c * k / (Hp * dtaudx) * Theta_p[l - 1];
    }
  }

  // // SET: Neutrino perturbations (N_ell)
  if (neutrinos)
  {
    Nu[0] = Nu_tc[0];
    Nu[1] = Nu_tc[1];
    Nu[2] = Nu_tc[2];
    for (int l = 3; l < n_ell_neutrinos; l++)
    {
      Nu[l] = Nu_tc[l];
    }
  }

  // check for nans
  for (int i = 0; i < y.size(); i++)
  {
    if (std::isnan(y[i]))
    {
      printf("y[%d] is nan\n", i);
      throw std::runtime_error("y is nan");
    }
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const
{
  double x_tight_coupling_end;

  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  for (int i = 0; i < n_x; i++)
  {
    double x = x_array[i];
    if (x > -8.3)
    {
      x_tight_coupling_end = -8.3;
      break;
    }
    double dtaudx = rec->dtaudx_of_x(x);
    double Hp = cosmo->Hp_of_x(x);
    if (abs(dtaudx) < 10. * std::min(1., c * k / Hp))
    {
      x_tight_coupling_end = x;
      break;
    }
  }

  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions()
{
  Utils::StartTiming("source");

  Vector k_array(n_k);
  k_array = Utils::logspace(log10(k_min), log10(k_max), n_k);
  Vector x_array(n_x);
  x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for (auto ix = 0; ix < x_array.size(); ix++)
  {
    const double x = x_array[ix];
    for (auto ik = 0; ik < k_array.size(); ik++)
    {
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      double eta_0 = cosmo->eta_of_x(0.0);
      double eta = cosmo->eta_of_x(x);
      double Hp = cosmo->Hp_of_x(x);
      double dHpdx = cosmo->dHpdx_of_x(x);
      double ddHpdxx = cosmo->ddHpddx_of_x(x);

      double g_tilde = rec->g_tilde_of_x(x);
      double dg_tildedx = rec->dgdx_tilde_of_x(x);
      double ddg_tildedxx = rec->ddgddx_tilde_of_x(x);
      double tau = rec->tau_of_x(x);

      double dPhidx = Phi_spline.deriv_x(x, k);
      double Psi = get_Psi(x, k);
      double dPsidx = Psi_spline.deriv_x(x, k);
      double Pi = get_Pi(x, k);
      double dPidx = Pi_spline.deriv_x(x, k);
      double ddPidxx = Pi_spline.deriv_xx(x, k);
      double v_b = get_v_b(x, k);
      double dv_bdx = v_b_spline.deriv_x(x, k);
      Vector Theta(Constants.n_ell_theta);
      for (int l = 0; l < Constants.n_ell_theta; l++)
      {
        Theta[l] = get_Theta(x, k, l);
      }

      // Temperatur source
      ST_array[index] = g_tilde * (Theta[0] + Psi + 1. / 4. * Pi) + exp(-tau) * (dPsidx - dPhidx) - 1. / (c * k) * (dHpdx * g_tilde * v_b + Hp * dg_tildedx * v_b + Hp * g_tilde * dv_bdx) + 3. / (4. * c * c * k * k) * (dHpdx * (dHpdx * g_tilde * Pi + Hp * dg_tildedx * Pi + Hp * g_tilde * dPidx) + Hp * (ddHpdxx * g_tilde * Pi + 2. * dHpdx * dg_tildedx * Pi + 2. * dHpdx * g_tilde * dPidx + Hp * ddg_tildedxx * Pi + 2. * Hp * dg_tildedx * dPidx + Hp * g_tilde * ddPidxx));

      // Polarization source
      if (Constants.polarization)
      {
        SE_array[index] = 3. * g_tilde * Pi / (4. * k * k * (eta_0 - eta) * (eta_0 - eta));
      }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  if (Constants.polarization)
  {
    SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx)
{

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;
  const bool neutrinos = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm = y[Constants.ind_deltacdm_tc];
  const double &delta_b = y[Constants.ind_deltab_tc];
  const double &v_cdm = y[Constants.ind_vcdm_tc];
  const double &v_b = y[Constants.ind_vb_tc];
  const double &Phi = y[Constants.ind_Phi_tc];
  const double *Theta = &y[Constants.ind_start_theta_tc];
  const double *Nu = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx = dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx = dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx = dydx[Constants.ind_vcdm_tc];
  double &dv_bdx = dydx[Constants.ind_vb_tc];
  double &dPhidx = dydx[Constants.ind_Phi_tc];
  double *dThetadx = &dydx[Constants.ind_start_theta_tc];
  double *dNudx = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  double a = exp(x);

  double Hp = cosmo->Hp_of_x(x);
  double dHpdx = cosmo->dHpdx_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);
  double ddtauddx = rec->ddtauddx_of_x(x);
  double eta = cosmo->eta_of_x(x);

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi = -Phi - 12. * H0 * H0 / (c * c * k * k * a * a) * Omega_nu0 * Nu[2];
  // check for nans
  if (std::isnan(Psi))
  {
    printf("Psi is nan in TC\n");
    throw std::runtime_error("Psi is nan");
  }
  double R = 4. * Omega_gamma0 / (3. * Omega_b0 * a);
  double q = (-((1. - R) * dtaudx + (1 + R) * ddtauddx) * (3. * Theta[1] + v_b) - c * k / Hp * Psi + (1. - dHpdx / Hp) * c * k / Hp * (-Theta[0]) - c * k / Hp * dThetadx[0]) / ((1. + R) * dtaudx + dHpdx / Hp - 1.);

  dPhidx = Psi - c * c * k * k / (3. * Hp * Hp) * Phi + H0 * H0 / (2. * Hp * Hp) * (Omega_CDM0 / a * delta_cdm + Omega_b0 / a * delta_b + 4. * Omega_gamma0 / (a * a) * Theta[0] + 4. * Omega_nu0 / (a * a) * Nu[0]);
  ddelta_cdmdx = c * k / Hp * v_cdm - 3. * dPhidx;
  dv_cdmdx = -v_cdm - c * k / Hp * Psi;
  ddelta_bdx = c * k / Hp * v_b - 3. * dPhidx;
  dv_bdx = 1. / (1. + R) * (-v_b - c * k / Hp * Psi + R * (q + c * k / Hp * (-Theta[0]) - c * k / Hp * Psi));

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -c * k / Hp * Theta[1] - dPhidx;
  dThetadx[1] = 1. / 3. * (q - dv_bdx);

  // SET: Neutrino mutlipoles (Nu_ell)
  if (neutrinos)
  {
    dNudx[0] = -c * k / Hp * Nu[1] - dPhidx;
    dNudx[1] = c * k / (3. * Hp) * Nu[0] - 2. * c * k / (3. * Hp) * Nu[2] + c * k / (3. * Hp) * Psi;
    for (int l = 2; l < n_ell_neutrinos_tc - 1; l++)
    {
      dNudx[l] = l * c * k / ((2. * l + 1.) * Hp) * Nu[l - 1] - (l + 1.) * c * k / ((2 * l + 1.) * Hp) * Nu[l + 1];
    }
    dNudx[n_ell_neutrinos_tc - 1] = c * k / Hp * Nu[n_ell_neutrinos_tc - 2] - c * (n_ell_neutrinos_tc) / (Hp * eta) * Nu[n_ell_neutrinos_tc - 1];
  }

  // check for nans
  if (std::isnan(dPhidx))
  {
    printf("dPhidx is nan\n");
    throw std::runtime_error("dPhidx is nan");
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx)
{

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta = Constants.n_ell_theta;
  const int n_ell_thetap = Constants.n_ell_thetap;
  const int n_ell_neutrinos = Constants.n_ell_neutrinos;
  const bool polarization = Constants.polarization;
  const bool neutrinos = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm = y[Constants.ind_deltacdm];
  const double &delta_b = y[Constants.ind_deltab];
  const double &v_cdm = y[Constants.ind_vcdm];
  const double &v_b = y[Constants.ind_vb];
  const double &Phi = y[Constants.ind_Phi];
  const double *Theta = &y[Constants.ind_start_theta];
  const double *Theta_p = &y[Constants.ind_start_thetap];
  const double *Nu = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx = dydx[Constants.ind_deltacdm];
  double &ddelta_bdx = dydx[Constants.ind_deltab];
  double &dv_cdmdx = dydx[Constants.ind_vcdm];
  double &dv_bdx = dydx[Constants.ind_vb];
  double &dPhidx = dydx[Constants.ind_Phi];
  double *dThetadx = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx = &dydx[Constants.ind_start_thetap];
  double *dNudx = &dydx[Constants.ind_start_nu];

  double a = exp(x);

  // Cosmological parameters and variables
  double Hp = cosmo->Hp_of_x(x);
  double eta = cosmo->eta_of_x(x);

  // Recombination variables
  double dtaudx = rec->dtaudx_of_x(x);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // check phi for nans
  if (std::isnan(Phi))
  {
    printf("Phi is nan\n");
    throw std::runtime_error("Phi is nan");
  }

  double Psi = -Phi - 12. * H0 * H0 / (c * c * k * k * a * a) * (Omega_gamma0 * Theta[2] + Omega_nu0 * Nu[2]);
  // check for nans
  if (std::isnan(Psi))
  {
    printf("Psi is nan\n");
    throw std::runtime_error("Psi is nan");
  }
  double R = 4. * Omega_gamma0 / (3. * Omega_b0 * a);
  double Pi = Theta[2] + Theta_p[0] + Theta_p[2];

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - c * c * k * k / (3. * Hp * Hp) * Phi + H0 * H0 / (2. * Hp * Hp) * (Omega_CDM0 / a * delta_cdm + Omega_b0 / a * delta_b + 4. * Omega_gamma0 / (a * a) * Theta[0] + 4. * Omega_nu0 / (a * a) * Nu[0]);
  ddelta_cdmdx = c * k / Hp * v_cdm - 3. * dPhidx;
  dv_cdmdx = -v_cdm - c * k / Hp * Psi;
  ddelta_bdx = c * k / Hp * v_b - 3. * dPhidx;
  dv_bdx = -v_b - c * k / Hp * Psi + dtaudx * R * (3 * Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -c * k / Hp * Theta[1] - dPhidx;
  dThetadx[1] = c * k / (3. * Hp) * Theta[0] - 2. * c * k / (3. * Hp) * Theta[2] + c * k / (3. * Hp) * Psi + dtaudx * (Theta[1] + 1. / 3. * v_b);
  for (int l = 2; l < n_ell_theta - 1; l++)
  {
    dThetadx[l] = l * c * k / ((2. * l + 1.) * Hp) * Theta[l - 1] - (l + 1) * c * k / ((2. * l + 1.) * Hp) * Theta[l + 1] + dtaudx * (Theta[l]);
    if (l == 2)
    {
      dThetadx[2] = dThetadx[2] - dtaudx * 1. / 10. * Pi;
    }
  }
  dThetadx[n_ell_theta - 1] = c * k / Hp * Theta[n_ell_theta - 2] - c * (n_ell_theta) / (Hp * eta) * Theta[n_ell_theta - 1] + dtaudx * Theta[n_ell_theta - 1];

  // SET: Photon polarization multipoles (Theta_p_ell)
  if (polarization)
  {
    dTheta_pdx[0] = -c * k / Hp * Theta_p[1] + dtaudx * (Theta_p[0] - 1. / 2. * Pi);
    for (int l = 1; l < n_ell_thetap - 1; l++)
    {
      dTheta_pdx[l] = l * c * k / ((2. * l + 1.) * Hp) * Theta_p[l - 1] - (l + 1) * c * k / ((2. * l + 1.) * Hp) * Theta_p[l + 1] + dtaudx * Theta_p[l];
      if (l == 2)
      {
        dTheta_pdx[2] = dTheta_pdx[2] - dtaudx * 1. / 10. * Pi;
      }
    }
    dTheta_pdx[n_ell_thetap - 1] = c * k / Hp * Theta_p[n_ell_thetap - 2] - c * (n_ell_thetap) / (Hp * eta) * Theta_p[n_ell_thetap - 1] + dtaudx * Theta_p[n_ell_thetap - 1];
  }

  // // SET: Neutrino mutlipoles (Nu_ell)
  if (neutrinos)
  {
    dNudx[0] = -c * k / Hp * Nu[1] - dPhidx;
    dNudx[1] = c * k / (3. * Hp) * Nu[0] - 2. * c * k / (3. * Hp) * Nu[2] + c * k / (3. * Hp) * Psi;
    for (int l = 2; l < n_ell_neutrinos - 1; l++)
    {
      dNudx[l] = l * c * k / ((2. * l + 1.) * Hp) * Nu[l - 1] - (l + 1.) * c * k / ((2. * l + 1.) * Hp) * Nu[l + 1];
    }
    dNudx[n_ell_neutrinos - 1] = c * k / Hp * Nu[n_ell_neutrinos - 2] - c * (n_ell_neutrinos) / (Hp * eta) * Nu[n_ell_neutrinos - 1];
  }

  // find which parameter is nan
  for (int i = 0; i < Constants.n_ell_tot_full; i++)
  {
    if (std::isnan(y[i]))
    {
      printf("y[%d] is nan\n", i);
      throw std::runtime_error("y is nan");
    }
  }

  // check for nans
  for (int i = 0; i < Constants.n_ell_tot_full; i++)
  {
    if (std::isnan(dydx[i]))
    {
      printf("dydx[%d] is nan\n", i);
      throw std::runtime_error("dydx is nan");
    }
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const
{
  return delta_cdm_spline(x, k);
}
double Perturbations::get_delta_b(const double x, const double k) const
{
  return delta_b_spline(x, k);
}
double Perturbations::get_v_cdm(const double x, const double k) const
{
  return v_cdm_spline(x, k);
}
double Perturbations::get_v_b(const double x, const double k) const
{
  return v_b_spline(x, k);
}
double Perturbations::get_Phi(const double x, const double k) const
{
  return Phi_spline(x, k);
}
double Perturbations::get_Psi(const double x, const double k) const
{
  return Psi_spline(x, k);
}
double Perturbations::get_Pi(const double x, const double k) const
{
  return Pi_spline(x, k);
}
double Perturbations::get_Source_T(const double x, const double k) const
{
  return ST_spline(x, k);
}
double Perturbations::get_Source_E(const double x, const double k) const
{
  return SE_spline(x, k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const
{
  return Theta_spline[ell](x, k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const
{
  return Theta_p_spline[ell](x, k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const
{
  return Nu_spline[ell](x, k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const
{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start << "\n";
  std::cout << "x_end:         " << x_end << "\n";
  std::cout << "n_x:     " << n_x << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc << "\n";
  std::cout << "n_k:     " << n_k << "\n";
  if (Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if (Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta << "\n";
  if (Constants.polarization)
  {
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap << "\n";
  }
  if (Constants.neutrinos)
  {
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc << "\n";
  if (Constants.neutrinos)
  {
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const
{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&](const double x)
  {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x << " ";
    fp << get_Theta(x, k, 0) << " ";
    fp << get_Theta(x, k, 1) << " ";
    fp << get_Theta(x, k, 2) << " ";
    fp << get_Phi(x, k) << " ";
    fp << get_Psi(x, k) << " ";
    fp << get_Pi(x, k) << " ";
    fp << get_Source_T(x, k) << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(5, arg) << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(50, arg) << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(500, arg) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

extern "C"
{
  Perturbations *Perturbations_new(BackgroundCosmology *cosmo, RecombinationHistory *rec)
  {
    return new Perturbations(cosmo, rec);
  }
  void Perturbations_solve(Perturbations *pert)
  {
    pert->solve();
  }
  double Perturbations_get_delta_cdm(Perturbations *pert, const double x, const double k)
  {
    return pert->get_delta_cdm(x, k);
  }
  double Perturbations_get_delta_b(Perturbations *pert, const double x, const double k)
  {
    return pert->get_delta_b(x, k);
  }
  double Perturbations_get_Theta(Perturbations *pert, const double x, const double k, const double ell)
  {
    return pert->get_Theta(x, k, ell);
  }
  double Perturbations_get_Nu(Perturbations *pert, const double x, const double k, const double ell)
  {
    return pert->get_Nu(x, k, ell);
  }
  double Perturbations_get_v_cdm(Perturbations *pert, const double x, const double k)
  {
    return pert->get_v_cdm(x, k);
  }
  double Perturbations_get_v_b(Perturbations *pert, const double x, const double k)
  {
    return pert->get_v_b(x, k);
  }
  double Perturbations_get_Phi(Perturbations *pert, const double x, const double k)
  {
    return pert->get_Phi(x, k);
  }
  double Perturbations_get_Psi(Perturbations *pert, const double x, const double k)
  {
    return pert->get_Psi(x, k);
  }
  double Perturbations_get_Theta_p(Perturbations *pert, const double x, const double k, const double ell)
  {
    return pert->get_Theta_p(x, k, ell);
  }
}