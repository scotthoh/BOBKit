.. highlight: python

Useful functions
================

BuccaneerLog class
------------------
.. autosummary::

   bobkit.buccaneer.Log

This class can help with logging and profiling of the various Buccaneer steps.
It consists of methods to:

* write xml output of the summary from each Buccaneer cycle
* profile the time taken for each of the step logged
* output a summary of residues built and sequenced as string

For example

.. doctest::

   >>> from bobkit.buccaneer import Log, Ca_join
   >>> log = Log("test")
   >>> cajoin = Ca_join(2.0, 2.0)
   >>> cajoin(mol_wrk)
   >>> # log the time for joining step
   >>> log.log("JOIN", mol_wrk, False)
   >>> # write xml output
   >>> log.xml("output.xml")
   >>> # output a summary of residues built and sequenced as string
   >>> summary = log.summary(mol_wrk,mol_mr, seq_wrk)
   >>> print(summary)
   126 residues were built in 126 fragments, the longest having    1 residues.
         0 residues were sequenced, after pruning.
         0 residues were uniquely allocated to   0 chains.
     Completeness by residues built:     0.0%
     Completeness of chains (number):    0.0%    (0)
   >>> # example output of profile
   >>> log.profile()
   Profile:
   JOIN:     0.58 s


Setting reference for Buccaneer
-------------------------------
.. autosummary::

   bobkit.buccaneer.set_reference

``bobkit.buccaneer.set_reference()`` can be used to set file paths to reference data (MTZ & PDB).

.. doctest::

   >>> from bobkit.buccaneer import set_reference
   >>> mtz_ref = ""
   >>> pdb_ref = ""
   >>> set_reference(mtz_ref, pdb_ref)


Utility submodule
-----------------
The utility submodule also consist of methods to read and write structure file(s).
If ``bobkit.util.read_structure()`` is called without passing it a 
``bobkit.clipper.MiniMol`` instance, it will return a MiniMol instance containing 
the structure read from the file given. Otherwise, the MiniMol instance is updated 
with the strcuture read.

.. doctest::

   >>> from bobkit.clipper import MiniMol
   >>> from bobkit.util import read_structure as read_structure
   >>> # Reading structure without passing in a MiniMol instance
   >>> mmol = read_structure("test_data/pdb5ni1_cryst1.pdb")
   PDB file: ../../test_data/pdb5ni1_cryst1.pdb
     Number of atoms read: 4579
        0     N    45.716   55.727   67.167
     4578     O    66.745   51.174   62.217
   >>> # Reading structure by passing it a MiniMol instance
   >>> mmol2 = MiniMol()
   >>> flag = Util.read_structure(mmol2, "test_data/pdb5ni1_cryst1.pdb")
   PDB file: ../../test_data/pdb5ni1_cryst1.pdb
     Number of atoms read: 4579
        0     N    45.716   55.727   67.167
     4578     O    66.745   51.174   62.217

To write a pdb/cif file from MiniMol, use ``bobkit.util.write_structure()``.
The corresponding file extension will be assigned within method.

.. doctest::
   
   >>> from bobkit.util import write_structure as write_structure
   >>> # writing a pdb file
   >>> write_structure(mmol, "write_out.pdb", cif_format=False)
   >>> # writing a cif file
   >>> write_structure(mmol, "write_out.cif", cif_format=True)