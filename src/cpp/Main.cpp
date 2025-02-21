#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv)
{
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h = 0.67;
  double OmegaB0 = 0.05;
  double OmegaCDM0 = 0.267;
  double OmegaK0 = 0.0;
  double Neff = 3.046;
  double TCMB = 2.7255;

  // Recombination parameters
  double Yp = 0.245;

  // Power-spectrum parameters
  double A_s = 2.1e-9;
  double n_s = 0.965;
  double kpivot_mpc = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB);
  cosmo.solve();
  cosmo.info();

  // Output background evolution quantities
  cosmo.output("output/data/cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  mcmc_fit_to_supernova_data("data/supernovadata.txt", "output/data/results_supernovafitting.txt");

  // Print the cosmic time for x = -8.657694814439557
  std::cout << "Cosmic time at x = -8.657694814439557: " << cosmo.get_cosmic_time(-8.657694814439557) << " seconds" << std::endl;

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module II
  //=========================================================================

  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();

  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");

  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
