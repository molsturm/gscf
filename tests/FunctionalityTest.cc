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
class PlainScf
      : public gscf::PlainScf<
              PlainScfState<FockType, typename FockType::stored_matrix_type>> {
public:
  typedef FockType fock_type;
  typedef typename FockType::stored_matrix_type matrix_type;
  typedef gscf::PlainScf<PlainScfState<FockType, matrix_type>> base_type;
  typedef typename matrix_type::vector_type vector_type;
  typedef typename matrix_type::scalar_type scalar_type;
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
class PulayDiisScf : public gscf::PulayDiisScf<PulayDiisScfState<
                           FockType, typename FockType::stored_matrix_type>> {
public:
  typedef FockType fock_type;
  typedef typename fock_type::stored_matrix_type matrix_type;
  typedef gscf::PulayDiisScf<PulayDiisScfState<FockType, matrix_type>>
        base_type;
  typedef typename base_type::state_type state_type;

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

  // The expected eigenvectors
  MultiVector<vector_type> evec_expected{
        {-1.1884746467844811, 0., -0.24288050447076578, 0., 0.,
         0.16813161525502301, 0., 0., 0., 0., 0., 0.016396171817405908, 0., 0.},
        {-1.064786764522789, 0., 0.8777407505081162, 0., 0.,
         -0.3081669311487548, 0., 0., 0., 0., 0., -0.028869768632584114, 0.,
         0.},
        {0., 0., 0., -3.69049544150639e-9, 0.8573394853277652, 0., 0., 0., 0.,
         0., 0., 0., 0.00002919059875836615, -0.6818863586007807},
        {0., 0.9857660367413854, 0., 0., 0., 0., 0., 0., 0.47777120131625944,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0.8573394853277649, 3.69049544151069e-9, 0., 0., 0., 0.,
         0., 0., 0., -0.6818863586007805, -0.000029190598758366127},
        {-0.5840485708369669, 0., 0.05174625401524502, 0., 0.,
         -1.0729001918355632, 0., 0., 0., 0., 0., -0.07137766077631158, 0., 0.},
        {0., 0., 0., 1.1728582480320243e-9, -0.27246685510597846, 0., 0., 0.,
         0., 0., 0., 0., 0.000045420745388700296, -1.0610192320620837},
        {0., -0.033706141181753996, 0., 0., 0., 0., 0., 0., 1.0949264340797673,
         0., 0., 0., 0., 0.},
        {0., 0., 0., -0.2724668551059781, -1.1728582480253274e-9, 0., 0., 0.,
         0., 0., 0., 0., -1.0610192320620837, -0.00004542074538870031},
        {0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.000043369360575274524,
         0.9999999990595501, 0., 0., 0.},
        {-0.0019206466502236202, 0., -0.011672197660675484, 0., 0.,
         0.06685683586559842, 0., 0., 0., 0., 0., -0.9976924548257627, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., -0.9999999990595487,
         0.00004336936057527441, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.}};

  // The expected energies:
  scalar_type exp_energy_1e_terms = -14.91401133311085;
  scalar_type exp_energy_coulomb = 5.268082025558712;
  scalar_type exp_energy_exchange = -1.736263121264314;
  scalar_type exp_energy_total = -11.38219242881645;

  // Since we only converge to 5e-7, allow errors up to this degree.
  const double basetol = 5e-7;

  SECTION("PlainSCF") {
    error_wrapped_solvers::PlainScf<decltype(fock)> scf;
    auto res = scf.solve(fock, idata.s_bb());

    CHECK(res.n_iter() < 15);

    // Check the eigenvalues
    const auto& evalues = *res.eigenvalues_ptr();
    for (size_t i = 0; i < eval_expected.size(); ++i) {
      CHECK(evalues[i] == numcomp(eval_expected[i]).tolerance(2. * basetol));
    }

    const auto& evectors = *res.eigenvectors_ptr();
    // TODO For comparing all of them one needs to take rotations
    //      inside degenerate subspaces into account
    for (size_t i = 0; i < n_alpha; ++i) {
      linalgwrap::adjust_phase(evectors[i], evec_expected[i]);
      CHECK(evectors[i] == numcomp(evec_expected[i]).tolerance(10. * basetol));
    }

    // Check the energies:
    auto hf_energies = res.problem_matrix_ptr()->energies();
    CHECK(hf_energies.energy_1e_terms ==
          numcomp(exp_energy_1e_terms).tolerance(basetol));
    CHECK(hf_energies.energy_coulomb ==
          numcomp(exp_energy_coulomb).tolerance(basetol));
    CHECK(hf_energies.energy_exchange ==
          numcomp(exp_energy_exchange).tolerance(basetol));
    CHECK(hf_energies.energy_total ==
          numcomp(exp_energy_total).tolerance(basetol));
  }  // PlainSCF

  SECTION("PulayDiisScf") {
    error_wrapped_solvers::PulayDiisScf<decltype(fock)> scf(
          {{"n_prev_steps", size_t(4)}});
    auto res = scf.solve(fock, idata.s_bb());

    CHECK(res.n_iter() < 8);

    // Check the eigenvalues
    const auto& evalues = *res.eigenvalues_ptr();
    for (size_t i = 0; i < eval_expected.size(); ++i) {
      CHECK(evalues[i] == numcomp(eval_expected[i]).tolerance(2. * basetol));
    }

    const auto& evectors = *res.eigenvectors_ptr();
    // TODO For comparing all of them one needs to take rotations
    //      inside degenerate subspaces into account
    for (size_t i = 0; i < n_alpha; ++i) {
      INFO("Comparing vector " + std::to_string(i));
      linalgwrap::adjust_phase(evectors[i], evec_expected[i]);
      CHECK(evectors[i] == numcomp(evec_expected[i]).tolerance(10. * basetol));
    }

    // Check the energies:
    auto hf_energies = res.problem_matrix_ptr()->energies();
    CHECK(hf_energies.energy_1e_terms ==
          numcomp(exp_energy_1e_terms).tolerance(basetol));
    CHECK(hf_energies.energy_coulomb ==
          numcomp(exp_energy_coulomb).tolerance(basetol));
    CHECK(hf_energies.energy_exchange ==
          numcomp(exp_energy_exchange).tolerance(basetol));
    CHECK(hf_energies.energy_total ==
          numcomp(exp_energy_total).tolerance(basetol));
  }  // PulayDiisScf

}  // TEST_CASE

}  // namespace test
}  // namespace gscf
