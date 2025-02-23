#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;
class BackgroundCosmology
{
private:
  // Cosmological parameters
  double h;         // Little h = H0/(100km/s/Mpc)
  double OmegaB0;   // Baryon density today
  double OmegaCDM0; // CDM density today
  double OmegaK0;   // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
  double Neff;      // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
  double TCMB;      // Temperature of the CMB today in Kelvin

  // Derived parameters
  double H0;           // The Hubble parameter today H0 = 100h km/s/Mpc
  double OmegaLambda0; // Dark energy density today
  double OmegaGamma0;  // Photon density today (follows from TCMB)
  double OmegaNu0;     // Neutrino density today (follows from TCMB and Neff)

  // Start and end of x-integration (can be changed)
  double x_start = Constants.x_start;
  double x_end = Constants.x_end;

  // Splines to be made
  Spline eta_of_x_spline{"eta"};
  Spline cosmic_time_spline{"cosmic time"};

public:
  // Constructors
  BackgroundCosmology() = delete;
  BackgroundCosmology(
      double h,
      double OmegaB0,
      double OmegaCDM0,
      double OmegaK0,
      double Neff,
      double TCMB);

  // Print some useful info about the class
  void info() const;

  // Do all the solving
  void solve();

  // Output some results to file
  void output(const std::string filename) const;

  // Get functions that we must implement
  double eta_of_x(double x) const;
  double H_of_x(double x) const;
  double Hp_of_x(double x) const;
  double dHpdx_of_x(double x) const;
  double ddHpddx_of_x(double x) const;
  double get_OmegaB(double x = 0.0) const;
  double get_OmegaGamma(double x = 0.0) const;
  double get_OmegaNu(double x = 0.0) const;
  double get_OmegaCDM(double x = 0.0) const;
  double get_OmegaLambda(double x = 0.0) const;
  double get_OmegaK(double x = 0.0) const;
  double get_H0() const;
  double get_h() const;
  double get_Neff() const;
  double get_TCMB(double x = 0.0) const;

  // Distance measures
  double get_luminosity_distance_of_x(double x) const;
  double get_comoving_distance_of_x(double x) const;
  double get_cosmic_time(double x) const;
  double get_z(double x) const;
};

extern "C"
{
  BackgroundCosmology *BackgroundCosmology_new(double h, double OmegaB0, double OmegaCDM0, double OmegaK0, double Neff, double TCMB);
  void BackgroundCosmology_solve(BackgroundCosmology *bc);
  double BackgroundCosmology_get_cosmic_time(BackgroundCosmology *bc, double x);
  double BackgroundCosmology_get_OmegaGamma(BackgroundCosmology *bc, double x);
  double BackgroundCosmology_get_OmegaNu(BackgroundCosmology *bc, double x);
  double BackgroundCosmology_get_H_of_x(BackgroundCosmology *bc, double x);
  double BackgroundCosmology_get_dHpdx_of_x(BackgroundCosmology *bc, double x);
  double BackgroundCosmology_eta_of_x(BackgroundCosmology *bc, double x);
}

#endif
