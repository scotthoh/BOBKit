.. highlight: python

Filtering fragments in poor density
===================================

In the filtering routine, the residues which have not been docked into the sequence
and are in poor density are removed. Fragments with sequences of less than 6 residues 
are removed.

Ca_filter class
---------------
.. autosummary::

   bobkit.buccaneer.Ca_filter
  
The constructor takes in an argument for the sigma cutoff (default: 3.0 Å, Buccaneer use 1.0 Å).
An instance of the class can be called by passing it the working model and a working map.

.. doctest::

  >>> from bobkit.buccaneer import Ca_filter
  >>> cafltr = Ca_filter(1.0)
  >>> cafltr(mol_wrk, xwrk)

There are also two static methods ``filter`` which can be used to run the filtering routine.
One takes in the following arguments:

* ``mol`` - working model
* ``xmap`` - working map
* ``sigcut`` - sigma cutoff to keep residues
* ``keep`` - boolean argument to keep or delete residues in poor density

.. doctest::

  >>> # keep=True will keep the residues in poor density
  >>> Ca_filter.filter(mol_wrk, xwrk, 1.0, keep=False)

or the other without the working map, in which the residues are scored based on their U-values.
Residues with good B factors will be kept. The arguments are:

* ``mol`` - working model
* ``sigcut`` - sigma cutoff to keep residues

  >>> Ca_filter.filter(mol_wrk, 1.0)

