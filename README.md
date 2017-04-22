# gscf

A library of algorithms for solving *non-linear generalised eigenproblems*
as they frequently occur when finding a *self-consistent field* (SCF) solution
to a differential equation in theoretical Physics or Chemistry.

At the moment the major application of the code is electronic structure theory,
more precisely solving the
[Hartree-Fock approximation](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method)
to the many-body Schrödinger equation.
Most available algorithms, however, are implemented in a way such that they are *not*
specific to this problem or the *type of basis* used to model the physical system.
Hence the name ``gscf`` for **g**eneralised **SCF**.

## Dependencies
``gscf`` depends on the following:
- [krims](https://linalgwrap.org/krims) for many basic utilities
  (GenMap, Exception handling, Subscription pointers, ...)
- [linalgwrap](https://linalgwrap.org/) as a flexible linear algebra backend.
  ``linalgwrap`` has further dependencies in order to perform the
  linear algebra.
  See their [README.md](https://github.com/linalgwrap/linalgwrap/README.md) for details.

Testing ``gscf`` further requires
- [Catch](https://github.com/philsquared/Catch/) for the testing environment
- [rapidcheck](https://github.com/emil-e/rapidcheck) for property-based testing

For building ``gscf`` (see [below](#building-gscf)) you only need the ``gscf``
source code and the build dependencies of linalgwrap, i.e. armadillo,
LAPACK and some BLAS library to be present on your system.
The other dependencies (including ``krims`` and ``linalgwrap``)
can be automatically downloaded during the build process
if you choose to do so (set ``AUTOCHECKOUT_MISSING_REPOS`` to ``ON``,
more below)

## Building ``gscf``
All compilers starting from ``clang-3.5`` and ``gcc-4.8`` should be able to build the code.
We strictly require ``C++11`` to build this library.

If you choose to build with the flag ``AUTOCHECKOUT_MISSING_REPOS`` set to ``ON``
all required dependencies will be automatically downloaded
and compiled alongside ``gscf``.

In order to build ``gscf`` with tests (recommended) run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON ..
cmake --build .
ctest
```

In order to build without tests run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON -DGSCF_ENABLE_TESTS=OFF -DLINALGWRAP_ENABLE_TESTS=OFF \
    -DKRIMS_ENABLE_TESTS=OFF ..
cmake --build .
```

## Short description of ``gscf``
In ``gscf`` we wish to solve a nonlinear SCF problem like
```
A(l,x) x = l x
```
where ``A`` is a matrix depending on its own eigenvalues ``l`` and
eigenvectors ``x`` implicitly.  

The following SCF algorithms are currently implemented for this task:
- [*Plain SCF*](src/gscf/PlainScf.hh):
  Start with an initial guess for ``A`` and diagonalise it to
  obtain hopefully improved set of eigenpairs ``(l,x)``.
  From these build a new matrix ``A`` and repeat.
- [*DIIS*](src/gscf/PulayDiisScf.hh):
  An algorithms based on the [DIIS algorithm](https://en.wikipedia.org/wiki/DIIS)
  by Peter Pulay.
- [*Truncated Optimal Damping*](src/gscf/TruncatedOptDampScf.hh):
  Based on the *Optimal damping algorithm* by Eric Cancès and Claude le Bris.
  This algorithm is specific to the Hartree-Fock problem only.

