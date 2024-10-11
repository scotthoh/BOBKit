.. highlight: python

Correcting sequence
===================

After sequence are assigned to chains, corrections of insertions and deletions
as identified during the sequencing step are performed. This is done by rebuilding
to add or delete a residue where possible.

Ca_correct class
----------------
.. autosummary::

   bobkit.buccaneer.Ca_correct

The :py:meth:`constructor <bobkit.buccaneer.Ca_correct.__init__>` takes in an argument
for torsion sampling (default: 12). An instance of the class can be 
:py:meth:`called <bobkit.buccaneer.Ca_correct.__call__>` by passing it the working model, 
working map, likelihood map targets for the amino acids and the known sequence.

.. doctest::

  >>> from bobkit.buccaneer import Ca_correct
  >>> cacor = Ca_correct(12)
  >>> cacor(molwrk, xwrk, llkcls, input_sequence)
  >>> # return number of corrected C-alphas
