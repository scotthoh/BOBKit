.. highlight: python

Pruning fragments
=================

In this stage, the fragments which provide inconsistent interpretations of the same electron density are examined.
Those clashing residues with lower scores, fragments with sequences less than 6 residues, or with worse density are removed.

Ca_prune class
--------------
.. autosummary::

   bobkit.buccaneer.Ca_prune
  
The constructor takes in an argument for radius cutoff (default: 3.0 Å).
An instance of the Ca_prune class can be called by passing it the following arguments:

* ``mol`` - working model
* ``xmap`` - working map

.. doctest::

  >>> from bobkit.buccaneer import Ca_prune
  >>> caprune = Ca_prune(3.0)
  >>> caprune(mol_wrk, xwrk)

The equivalent static method ``Ca_prune.prune()`` takes in three arguments:

* ``mol`` - working model
* ``xmap`` - working map
* ``rad`` - radius cutoff, default is 3.0 Å

.. doctest::
  
  >>> Ca_prune.prune(mol_wrk, xwrk, 3.0)

