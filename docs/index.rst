.. Hazel documentation master file, created by
   sphinx-quickstart on Tue Apr  5 19:26:23 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Hazel
=================================
Hazel (an acronym for HAnle and ZEeman Light) is a computer program for the 
synthesis and inversion of Stokes profiles caused by the joint action of atomic 
level polarization and the Hanle and Zeeman effects. It is based on the quantum 
theory of spectral line polarization, which takes into account rigorously all the 
relevant physical mechanisms and ingredients: optical pumping, atomic level 
polarization, level crossings and repulsions, Zeeman, Paschen-Back and Hanle effects. 
The code is written in standard Fortran 90. Its parameters are passed using four 
configuration files that can be manually edited. These configuration files are heavily 
commented, so that their edition should be an easy task. In any case, two front-ends 
coded in IDL are given as a part of the distribution in order to facilitate a user-friendly 
execution of the program. A parallel version of the code using Message Passing Interface 
(MPI) is also available as well as a Python wrapper.


Indices and tables
==================

* :doc:`installation`
* :doc:`manual`
* :doc:`equations`
* :doc:`ambiguities`
* :ref:`modindex`
* :ref:`search`
* :doc:`disclaimer`

.. toctree::
   :hidden:
   :maxdepth: 2

   installation
   inputFiles
   atomicModels
   inputOutput
   graphical
   equations
   ambiguities
   disclaimer

