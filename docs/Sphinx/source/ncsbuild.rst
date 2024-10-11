.. highlight: python

Building noncrystallographic (NCS) symmetry
===========================================

Any noncrystallographic symmetry (NCS) relationship found in the model are used
to extend existing chains by combining all of the NCS-related chains.

Ca_ncsbuild class
-----------------
.. autosummary::

   bobkit.buccaneer.Ca_ncsbuild

The :py:meth:`constructor <bobkit.buccaneer.Ca_ncsbuild.__init__>` takes in three arguments:
 
* ``reliability`` - sequence reliability
* ``rmsd`` - root mean square deviation
* ``nmin`` - minimum number of matches

An instance of the class can be :py:meth:`called <bobkit.buccaneer.Ca_ncsbuild.__call__>` by 
passing it the following arguments:

* ``mol`` - working model
* ``xmap`` - working map
* ``llktarget`` - log likelihood map targets for the amino acids
* ``seq`` - known sequence

.. doctest::

  >>> from bobkit.buccaneer import Ca_ncsbuild
  >>> cancsbuild = Ca_ncsbuild(0.95, 1.9, 12)
  >>> cancsbuild(mol_wrk, xwrk, llkcls, seq)

