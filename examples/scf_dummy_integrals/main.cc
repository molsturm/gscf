#include "FockMatrix.hh"
#include "integrals/IntegralsSturmian14.hh"
#include <gscf/PlainScf.hh>
//#include <gscf/PulayDiisScfSolver.hh>
#include <gscf/version.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>

namespace scf_dummy {
using namespace linalgwrap;
using namespace gscf;

ArmadilloMatrix<double> hack_guess(const ArmadilloMatrix<double>& s_bb) {
    // apply Löwdin normalisation to the basis functions
    //   - Diagonalise the overlap
    //   - Take 1/\sqrt{evals} at the diagonal
    //   - results in orthonormalised basis functions

    const arma::mat& s_bb_data = s_bb.data();

    // Diagonalize s_bb
    arma::vec eval;
    arma::mat evec;
    arma::eig_sym(eval, evec, s_bb_data);

    // take 1/sqrt( . ) for each eigenvalue:
    std::transform(std::begin(eval), std::end(eval), std::begin(eval),
                   [](double elem) { return 1. / std::sqrt(elem); });

    // construct new Armadillo matrix of guess vectors:
    arma::mat guess = evec * arma::diagmat(eval);
    // If number of orbitals is different from number of basis functions,
    // guess should be rectangular.
    //      int b = s_bb.n_rows;
    //      guess = slice(evec * diagmat(eval), b, 5,0,0);

    // Return it properly enwrapped:
    // Note, that our matrices are row-major, but since armadillo
    // is column-major, we need to transpose it first
    return ArmadilloMatrix<double>(guess.t());
}

/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param Z       Number of nuclei
 * \param k_exp   Sturmian exponent
 * \param n_alpha Number of alpha electrons
 * \param n_beta  Number of beta electrons
 */
void run_sturmian14(double Z, double k_exp, size_t n_alpha, size_t n_beta) {
    typedef IntegralsSturmian14 idata_type;
    typedef idata_type::matrix_type matrix_type;
    typedef idata_type::scalar_type scalar_type;
    typedef FockMatrix<idata_type> fock_type;

    // Define integral data:
    idata_type idata(Z, k_exp);

    // Use a plain zero-guess for now:
    matrix_type guess = hack_guess(idata.s_bb());

    // XXX debug
    genout.write(guess, "guess");
    genout.write(idata.s_bb(), "sbb");
    // XXX debug

    fock_type fock{n_alpha, n_beta, idata, guess};

    // Define SCF types:
    typedef PlainScf<fock_type> scf_type;
    typedef typename scf_type::scf_state_type scf_state_type;
    typedef typename scf_type::scf_control_type scf_control_type;

    // Allocate SCF objects:
    scf_control_type scf_control;
    scf_type scf(scf_control);

    scf.update_eigenpair_handler([&](const scf_state_type& s) {
        // Here we have the new eigenpairs and the old fock matrix.
        genout.write(*s.eigenvalues_ptr,
                     "evals" + std::to_string(s.n_iter_count()));
        std::cout << "   New orbital eigenvalues: " << std::endl;

        // Print orbital evals:
        std::ostream_iterator<scalar_type> out_it(std::cout, " ");
        std::cout << "        ";
        std::copy(s.eigenvalues_ptr->begin(), s.eigenvalues_ptr->end(), out_it);
        std::cout << std::endl;
    });

    scf.update_problem_matrix_handler([](const scf_state_type& s) {
        typedef fock_type::energies_type energies_type;
        energies_type hf_energies = s.problem_matrix_ptr->get_energies();

        std::cout << "   Current HF energies:" << std::endl
                  << "     E_kin    = " << hf_energies.energy_kinetic
                  << std::endl
                  << "     E_v0     = " << hf_energies.energy_elec_nuc_attr
                  << std::endl
                  << "     E_coul   = " << hf_energies.energy_coulomb
                  << std::endl
                  << "     E_xchge  = " << hf_energies.energy_exchange
                  << std::endl
                  << std::endl
                  << "     E_1e     = " << hf_energies.energy_1e_terms
                  << std::endl
                  << "     E_2e     = " << hf_energies.energy_2e_terms
                  << std::endl
                  << std::endl
                  << "     E_total  = " << hf_energies.energy_total
                  << std::endl;
    });

    scf.pre_step_handler([](const scf_state_type& s) {
        std::cout << std::endl
                  << "Starting SCF iteration " << s.n_iter_count() << std::endl;
    });

    scf.solve_assert(fock, idata.s_bb());
}
}  // namespace scf_dummy

int main() {
    std::cout << "gscf version: " << gscf::version::version_string()
              << std::endl
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
