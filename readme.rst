==========================================
OpenMC Monte Carlo Particle Transport Code
==========================================

The OpenMC project aims to provide a fully-featured Monte Carlo particle
transport code based on modern methods. It is a constructive solid geometry,
continuous-energy transport code that uses ACE format cross sections. The
project started under the Computational Reactor Physics Group at MIT.

Complete documentation on the usage of OpenMC is hosted on GitHub at
http://mit-crpg.github.io/openmc/. If you are interested in the project or would
like to help and contribute, please send a message to the OpenMC User's Group
`mailing list`_.

------------
Installation
------------

Detailed `installation instructions`_ can be found in the User's Guide.

---------------------------
XEON PHI (MIC) INSTRUCTIONS
---------------------------

Please follow these steps:
  1) build the host binary normally by setting "MIC = no" in the Makefile
  2) remove all the host object and .mod files by typing "make hostclean"
  3) set "MIC = yes" in the makefile and type "make" again.  

This creates two binaries: "openmc" will run on the host and "openmc.mic"
will run on the Xeon Phi device.

Set the following environment variables,
  export I_MPI_MIC=enable
  export I_MPI_MIC_POSTFIX=.mic
  export I_MPI_ROOT=/soft/compilers/intel/impi/4.1.3.048
  export I_MPI_FABRICS=shm:tcp  # or shm:dapl on InfiniBand

The application can then be run on the host and multiple MIC devices:
  mpirun -n 1 -host localhost path/to/src/openmc : \
         -n 1 -host mic0 path/to/src/openmc : \
         -n 1 -host mic1 path/to/src/openmc

---------------
Troubleshooting
---------------

If you run into problems compiling, installing, or running OpenMC, first check
the `Troubleshooting section`_ in the User's Guide. If you are not able to find
a solution to your problem there, please send a message to the User's Group
`mailing list`_.

--------------
Reporting Bugs
--------------

OpenMC is hosted on GitHub and all bugs are reported and tracked through the
Issues_ feature on GitHub. However, GitHub Issues should not be used for common
troubleshooting purposes. If you are having trouble installing the code or
getting your model to run properly, you should first send a message to the
User's Group `mailing list`_. If it turns out your issue really is a bug in the
code, an issue will then be created on GitHub. If you want to request that a
feature be added to the code, you may create an Issue on github.

-------
License
-------

OpenMC is distributed under the MIT/X license_.

.. _mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-users
.. _installation instructions: http://mit-crpg.github.io/openmc/usersguide/install.html
.. _Troubleshooting section: http://mit-crpg.github.io/openmc/usersguide/troubleshoot.html
.. _Issues: https://github.com/mit-crpg/openmc/issues
.. _license: http://mit-crpg.github.io/openmc/license.html
