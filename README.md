# AdaptiveSampler - Importance sampling for Monte Carlo Integration

## Author

* Richard Jones, University of Connecticut, Storrs, CT

## Description

AdaptiveSampler is a c++ class that implements a very general adaptive
algorithm for efficient computation of a Monte Carlo integral over a
multi-dimensional domain. Computation of Monte Carlo integrals on a
multi-dimensional domain are often beset by the problem that the support
for the integral comes from a very small volume of the total domain of
integration. In such a situation, the Monte Carlo estimation can be very
slow to converge. This class attempts to zero in on the "hot spots" in
the integration domain and reduce the variance on the result for a given
sample size by using importance sampling. The advantage of this class
over implementing your own is that the procedure it uses is strictly
empirical and does not require that the user provide any hints about
the behavior of the integrand. It simply requires that the domain of
integration be transformed by variable substitution into the unit 
hypercube (0,1)^N, and that the integral be convergent.

## History

This package was originally written to provide the compute the polarization
of coherent bremsstrahlung photon beam by integrating over the recoil
electron four-momentum. It was fairly simple to generalize into a tool for
any Monte Carlo integration problem which can be written as an integral
over the unit hypercube. 

## Release history

See VERSIONS file in the project directory.

## Usage synopsis

General functions that support plotting of cross sections and particle
generation are provided for specific QED processes in the source files
with names that end with ".C". These functions call methods of the
TCrossSection class to compute differential cross sections and rates.
See the comments within those source files for details.

## Dependencies

You need a working c++ compiler with support for c++11 language features.
To run the example program out of the box, you also need to install the
ROOT package from CERN. If you do not have ROOT already running on your
system, it will be simpler to modify the ex.cc function to replace the
TRandom object with a uniform(0,1) random number generator of your own.

## Building instructions

The following command compiles the example program ex.cc and builds an
executable ex.

    $ make

You can run ex without any command line arguments, or you can simply use
it as an example for how to build an application to do your Monte Carlo
integral.

## Documentation

See comments at the head of the following source files:

1. AdaptiveSampler.hh - contains usage explanations and an example.
2. AdaptiveSampler.cc - contains implementation notes by the author
                        and hints for how to use it effectively.
3. ex.cc - an example code with minimal explanation, useful as as
           starting point for your own projects.

## Troubleshooting

The code is open-source, and the behavior should be apparent from
the names of the methods. For help with bugs and questions about
how to implement additional processes, please contact the author.

## Bugs

None known at this time.

## How to contribute

## Contact the author

Write to richard.t.jones at uconn.edu.
