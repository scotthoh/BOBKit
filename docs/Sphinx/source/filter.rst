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
  
The :py:meth:`constructor <bobkit.buccaneer.Ca_filter.__init__>` takes in an argument 
for the sigma cutoff (default: 3.0 Å, Buccaneer use 1.0 Å). An instance of the class 
can be :py:meth:`called <bobkit.buccaneer.Ca_filter.__call__>` by passing it the working 
model and a working map.

.. doctest::

  >>> from bobkit.buccaneer import Ca_filter
  >>> cafltr = Ca_filter(1.0)
  >>> cafltr(mol_wrk, xwrk)

There is also a static method :py:meth:`filter <bobkit.buccaneer.Ca_filter.filter>` 
which can be used to run the filtering routine, with or without a working map.
If the method is called without giving it a working map, the residues are scored
based on their U-values. Residues with good B factors will be kept.

.. doctest::

  >>> # with working map, xwrk
  >>> # keep=True will keep the residues in poor density
  >>> Ca_filter.filter(mol_wrk, xwrk, 1.0, keep=False)
  >>>
  >>> # without working map
  >>> Ca_filter.filter(mol_wrk, 1.0)

