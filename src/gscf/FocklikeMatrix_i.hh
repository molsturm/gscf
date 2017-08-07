//
// Copyright (C) 2017 by the gscf authors
//
// This file is part of gscf.
//
// gscf is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gscf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gscf. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "ScfProblemMatrix_i.hh"
#include <krims/Range.hh>

namespace gscf {

/** Enum object which specifies the orbital space to select
 *  in a Fock-like SCF problem */
enum class OrbitalSpace {
  /** The occupied alpha orbitals */
  OCC_ALPHA,

  /** The occupied beta orbitals */
  OCC_BETA,

  /** The virtual (empty) alpha orbitals */
  VIRT_ALPHA,

  /** The virtual (empty) beta orbitals */
  VIRT_BETA,
};

/** This class defines the extra interface, which some SCF solvers
 *  require from the Problem matrix to be solved.
 *
 *  It essentially is an extension to the default, generic, problem
 *  matrix interface given by a lazyten::LazyMatrixExpression,
 *  which is required for some SCF solvers which are used especially
 *  in quantum chemistry.
 *
 *  TODO Keep this interface as small as possible. Replace by more
 *       general concepts if they are available.
 */
template <typename Scalar>
class FocklikeMatrix_i : public ScfProblemMatrix_i {
 public:
  typedef Scalar scalar_type;

  // TODO Some of the things here might go once we have the block diagonal
  //      matrix data structures

  /** Return the one electron energy, i.e. the energy of the
   *  terms independent of the electron density.
   *
   *  This energy is usually proportional to \f$tr(h D)\f$,
   *  where \f$h\f$ is the part of this matrix, which is
   *  constant wrt. changes in the density \f$D\f$ (the core
   *  hamiltonian).
   */
  virtual scalar_type energy_1e_terms() const = 0;

  /** Return the two electron energy.
   *
   *  This energy is usually proportional to \f$\frac 1 2 tr(G D)\f$,
   *  where \f$G\f$ is the part of this matrix, which is
   *  dependant on the density \f$D\f$, i.e. the Coulomb and Exchange
   *  matrices.
   */
  virtual scalar_type energy_2e_terms() const = 0;

  // TODO Conceptionally it would be better to return some sort of view
  //      mask here, which allows a more transparent handling of the
  //      alphas and betas

  /** Return the index range which is used by orbitals of the
   *  specified subspace.
   *
   * E.g. if oss == OCC_ALPHA the range will run over the indices of
   * all occupied alpha MOs. For restricted calculations it can be assumed
   * that OCC_ALPHA == OCC_BETA and VIRT_ALPHA == VIRT_BETA
   *
   * \note These ranges are only valid for coefficients and other objects
   * which have the same structure as the coefficients which were used
   * in the most recent call to update the Fock matrix. In other words
   * usually this applies only to quantities from the scf iteration in
   * which this very fock matrix was updated, i.e. the orbital energies
   * and the coeffients, which went into making this Fock matrix.
   * If the ranges are applied to other objects, specifically if
   * these other objects have a different number solution orbitals,
   * then the index ranges from here will not be applicable.
   */
  virtual krims::Range<size_t> indices_orbspace(OrbitalSpace osp) const = 0;

  //** Return a reference to the coefficients which were used to build
  // *  the current problem matrix.
  // *
  // *  These are exactly the coefficients, which have been supplied
  // *  by the most recent call to update() to this matrix object.
  // */
  // virtual const lazyten::MultiVector<vector_type>& coefficients() const = 0;

  FocklikeMatrix_i()                        = default;
  virtual ~FocklikeMatrix_i()               = default;
  FocklikeMatrix_i(const FocklikeMatrix_i&) = default;
  FocklikeMatrix_i(FocklikeMatrix_i&&)      = default;
  FocklikeMatrix_i& operator=(FocklikeMatrix_i&&) = default;
  FocklikeMatrix_i& operator=(const FocklikeMatrix_i&) = default;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a focklike matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef vector_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsFocklikeMatrix : public std::false_type {};

template <typename T>
struct IsFocklikeMatrix<T, krims::VoidType<typename T::scalar_type>>
      : public std::is_base_of<FocklikeMatrix_i<typename T::scalar_type>, T> {};
//@}

}  // namespace gscf
