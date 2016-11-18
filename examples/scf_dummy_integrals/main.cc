#include "FockMatrix.hh"
#include "ScfLibrary.hh"
#include "integrals/IntegralsSturmian14.hh"
#include "loewdin_guess.hh"
#include <gscf/version.hh>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>

namespace scf_dummy {
using namespace linalgwrap;
using namespace gscf;

/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param Z       Number of nuclei
 * \param k_exp   Sturmian exponent
 * \param n_alpha Number of alpha electrons
 * \param n_beta  Number of beta electrons
 */
void run_sturmian14(double Z, double k_exp, size_t n_alpha, size_t n_beta) {
  // Define types for fock operator
  typedef IntegralsSturmian14 idata_type;
  typedef idata_type::scalar_type scalar_type;
  typedef FockMatrix<idata_type> fock_type;

  // Output file for mathematica debug:
  std::ofstream mathematicafile("/tmp/debug_gscf_scfdummy_sturm14.m");
  auto debugout = linalgwrap::io::make_formatted_stream_writer<
        linalgwrap::io::Mathematica, scalar_type>(mathematicafile, 1e-5);

  // Define integral data:
  idata_type idata(Z, k_exp);

  // Write sbb to debug:
  debugout.write("sbb", idata.s_bb());

  // Use a Löwdin orthogonalised guess for now:
  auto guess = loewdin_guess(idata.s_bb());
  debugout.write("guess", guess);

  bool store_terms = true;
  fock_type fock{n_alpha, n_beta, idata, guess, store_terms};
  debugout.write("guessfock", fock);

  // Allocate SCF objects:
  // PlainScfHartreeFock<decltype(fock)> scfhf(debugout);
  DiisScfHartreeFock<decltype(fock)> scfhf(debugout);
  scfhf.update_control_params({{"n_prev_steps", size_t(4)}});
  scfhf.solve(fock, idata.s_bb());
}
}  // namespace scf_dummy

int main() {
  std::cout << "gscf version: " << gscf::version::version_string() << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl
            << std::endl;

  // This should be a test program using dummy integral data
  // and a simple fock matrix class in order to do some types of DIIS
  // and a plain SCF

  //
  // closed-shell Beryllium with k_exp = 1.351
  //
  double Z = 4.;
  double k_exp = 1.;  // 1.351;
  size_t n_alpha = 2;
  size_t n_beta = 2;
  scf_dummy::run_sturmian14(Z, k_exp, n_alpha, n_beta);

  return 0;
}
