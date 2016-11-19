#include <catch.hpp>
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <gscfmock/FockMatrix.hh>
#include <gscfmock/Integrals.hh>
#include <gscfmock/pulay_error.hh>
#include <linalgwrap/TestingUtils.hh>

namespace gscf {
namespace tests {
using namespace gscf;
using namespace linalgwrap;
using namespace krims;
using namespace gscfmock;

namespace error_wrapped_solvers {

template <typename FockType>
class PlainScf : public gscf::PlainScf<FockType, PlainScfState<FockType>> {
public:
  typedef FockType fock_type;
  typedef gscf::PlainScf<FockType, PlainScfState<FockType>> base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename matrix_type::vector_type vector_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::state_type state_type;

  //! The maximal allowed error for convergence
  scalar_type max_error = 5e-7;

  //! We converge if the frobenius norm of the Pulay error is
  // below above value
  bool is_converged(const state_type& s) const override {
    if (s.problem_matrix_ptr() == nullptr || s.eigenvectors_ptr() == nullptr)
      return false;

    const fock_type& fock_bb = *s.problem_matrix_ptr();
    const auto& coefficients_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();

    matrix_type error = pulay_error(fock_bb, coefficients_bf, overlap_bb);
    return norm_frobenius(error) < max_error;
  }

  PlainScf() : base_type() {}
  PlainScf(const ParameterMap& map) : base_type(map) {}
};

template <typename FockType>
class PulayDiisScf
      : public gscf::PulayDiisScf<FockType, PulayDiisScfState<FockType>> {
public:
  typedef FockType fock_type;
  typedef gscf::PulayDiisScf<FockType, PulayDiisScfState<FockType>> base_type;
  typedef typename base_type::state_type state_type;
  typedef typename state_type::matrix_type matrix_type;

  //! Define how we compute the error: Pulay error
  matrix_type calculate_error(const state_type& s) const override {
    const fock_type& fock_bb = *s.problem_matrix_ptr();
    const auto& coefficients_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();
    return pulay_error(fock_bb, coefficients_bf, overlap_bb);
  }

  PulayDiisScf() : base_type() {}
  PulayDiisScf(const ParameterMap& map) : base_type(map) {}
};

}  // namespace error_wrapped_solvers

TEST_CASE("SCF functionality test", "[SCF functionality]") {
  // The test problem
  double Z = 4.;  // Be atom
  double k_exp = 1.;
  size_t n_alpha = 2;
  size_t n_beta = 2;

  // Setup integral data
  typedef IntegralsSturmian14 idata_type;
  typedef typename idata_type::matrix_type matrix_type;
  typedef typename idata_type::scalar_type scalar_type;
  typedef typename matrix_type::vector_type vector_type;
  idata_type idata(Z, k_exp);

  // The guess to use ... we could use something
  // random here
  MultiVector<vector_type> guess(idata.nbas(), 4);
  guess[0](0) = -0.9238795325112872;
  guess[0](1) = -1.306562964876376;
  guess[0](5) = -0.9238795325112864;
  guess[1](3) = guess[1](7) = -0.9192110607898044;
  guess[2](2) = guess[2](6) = -0.9192110607898044;
  guess[3](4) = guess[3](8) = -0.9192110607898044;

  // The initial fock matrix
  FockMatrix<decltype(idata)> fock(n_alpha, n_beta, idata, guess);

  // The expected eigenvalues
  std::vector<scalar_type> eval_expected{
        -3.749602941400554,  -0.1755834257634913, -0.03942374133684358,
        0.08845442862852368, 0.08845442862852368, 0.2470759037725603,
        0.3644267886707621,  0.3644267886707621,  0.366327645695445,
        0.3727105903195644,  0.3727105903195644,  0.3761196059870835,
        0.3835526231630215,  0.3835526231630215};

  // The expected energies:
  scalar_type exp_energy_1e_terms = -14.91401133311085;
  scalar_type exp_energy_coulomb = 5.268082025558712;
  scalar_type exp_energy_exchange = -1.736263121264314;
  scalar_type exp_energy_total = -11.38219242881645;

  // Since we only converge to 5e-7, allow larger errors:
  auto highertol = NumCompConstants::change_temporary(
        1e-7 / std::numeric_limits<scalar_type>::epsilon());

  SECTION("PlainSCF") {
    error_wrapped_solvers::PlainScf<decltype(fock)> scf;
    auto res = scf.solve(fock, idata.s_bb());

    CHECK(res.n_iter_count() < 15);

    // Check the eigenvalues
    const auto evaltol = NumCompAccuracyLevel::Lower;
    const auto& evalues = *res.eigenvalues_ptr();
    for (size_t i = 0; i < eval_expected.size(); ++i) {
      CHECK(evalues[i] == numcomp(eval_expected[i]).tolerance(evaltol));
    }

    // Check the energies:
    const auto energytol = NumCompAccuracyLevel::Default;
    auto hf_energies = res.problem_matrix_ptr()->energies();
    CHECK(hf_energies.energy_1e_terms ==
          numcomp(exp_energy_1e_terms).tolerance(energytol));
    CHECK(hf_energies.energy_coulomb ==
          numcomp(exp_energy_coulomb).tolerance(energytol));
    CHECK(hf_energies.energy_exchange ==
          numcomp(exp_energy_exchange).tolerance(energytol));
    CHECK(hf_energies.energy_total ==
          numcomp(exp_energy_total).tolerance(energytol));
  }  // PlainSCF

  SECTION("PulayDiisScf") {
    error_wrapped_solvers::PulayDiisScf<decltype(fock)> scf(
          {{"n_prev_steps", size_t(4)}});
    auto res = scf.solve(fock, idata.s_bb());

    CHECK(res.n_iter_count() < 10);

    // Check the eigenvalues
    const auto evaltol = NumCompAccuracyLevel::Lower;
    const auto& evalues = *res.eigenvalues_ptr();
    for (size_t i = 0; i < eval_expected.size(); ++i) {
      CHECK(evalues[i] == numcomp(eval_expected[i]).tolerance(evaltol));
    }

    // Check the energies:
    const auto energytol = NumCompAccuracyLevel::Default;
    auto hf_energies = res.problem_matrix_ptr()->energies();
    CHECK(hf_energies.energy_1e_terms ==
          numcomp(exp_energy_1e_terms).tolerance(energytol));
    CHECK(hf_energies.energy_coulomb ==
          numcomp(exp_energy_coulomb).tolerance(energytol));
    CHECK(hf_energies.energy_exchange ==
          numcomp(exp_energy_exchange).tolerance(energytol));
    CHECK(hf_energies.energy_total ==
          numcomp(exp_energy_total).tolerance(energytol));
  }  // PulayDiisScf

}  // TEST_CASE

}  // namespace test
}  // namespace gscf
