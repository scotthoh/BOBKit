.. BOBKit documentation master file, created on Tue Sep 10 15:05:18 2024.

BOBKit - Model Building Kit in Python
=====================================
**B**\ iomacromolecular M\ **o**\ del **B**\ uilding **K**\ it

The BOBKit module provides python bindings the protein model building
software, Buccaneer_. Since some input arguments of Clipper types are 
needed for the Buccaneer_ methods, python bindings to some parts of the 
Clipper C++ libraries for X-ray crystallography computation are also
available in this module. This python module can aid developers who are 
developing machine learning softwares for model building.

The calculation in Buccaneer_ involves 10 stages, namely:

:doc:`Finding C-alphas <find>`
   Candidate C-alpha positions are located by searching the electron density.

:doc:`Growing chain fragments <grow>`
   The candidate C-alphas or input chains are grown by adding residues at
   either end, according to the density.

:doc:`Joining chain fragments <join>`
   Overlapping fragments are joined to make longer chains. If this leads
   to a junction in a chain, the contested residue is removed.

:doc:`Linking chain fragments <link>`
   Nearby N and C termini are examined to see if they can be linked by a
   short loop.

:doc:`Assigning sequence <sequence>`
   Likelihood comparison between the dnesity of each residue in the work
   structure and the residues of the reference structure allows sequence
   to be assigned to longer formats.

:doc:`Correcting sequence <correct>`
   Insertions and deletions in the model building are fixed by rebuilding,
   where possible.

:doc:`Filtering fragments in poor density <filter>`
   Residues in poor density are removed.

:doc:`Building NCS <ncsbuild>`
   Any NCS relationships found in the model are used to augment the
   related chains.

:doc:`Pruning fragments <prune>`
   Clashing fragments are examined and the one with the worse density is 
   removed.

:doc:`Rebuilding <rebuild>`
   Allows side chain atoms and carbonyl oxygens to be rebuilt.

These stages are exposed in python as classes. The associated methods
can be reused individually in various stages of any model building
workflows, provided that they given the correct input arguments.

BOBKit is open source and portable, it has been tested on Linux (Ubuntu 20.04).
The python bindings are written for Python 3.

This project is funded by the BBSRC Grant BB/X006492/1.

.. _Buccaneer: https://www.ccp4.ac.uk/html/cbuccaneer.html

.. note::

   This project is under development.
   Code refactoring might happen in newer versions.

Contents
--------

.. toctree::
   :caption: Getting started
   :maxdepth: 2
   
   What is BOBKit? <self>
   install
   pybuccaneer

.. toctree::
   :caption: What can BOBKit do?
   :maxdepth: 2

   Preparing log-likelihood map targets <prep>
   buccaneer_util

.. toctree::
   :caption: The 10 stages of Buccaneer
   :maxdepth: 2

   find
   grow
   join
   link
   sequence
   correct
   filter
   Building NCS <ncsbuild>
   prune
   Rebuilding <rebuild>

.. toctree::
   :caption: API Reference
   :maxdepth: 5

   py_api_reference
