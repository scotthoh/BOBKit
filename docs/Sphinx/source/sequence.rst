.. highlight: python

Assigning sequence
==================

In this routine, the likelihood map target is used to score each residue for being a particular type.
This is done by comparing the density of each residue in the working structure with those from of a reference
structure. Then, sequence alignment is performed between the given chain scores and the given sequence.
This allows longer fragments to be matched to the sequence.

Class which can be used for the sequencing routine:

* ``Ca_sequence`` for sequencing C-alpha chains using density
* ``Sequence_score_threaded`` for running multi-threaded sequence scoring routine
* ``Sequence_threaded`` for running multi-threaded sequencing routine

Ca_sequence class
-----------------
.. autosummary::

   bobkit.buccaneer.Ca_sequence

The constructor takes in an argument for sequence reliability (default: 0.5, Buccaneer use 0.95).

* ``Ca_sequence.set_cpus()`` can be used to set the number of cpu threads to use in running the calculations.
* ``Ca_sequence.set_semet()`` is used to set the flag to translate methionine (MET) to seleno-methionine (MSE).
* ``Ca_sequence.seq_prior_model()`` can be used to set a prior model for sequencing.

.. doctest::

  >>> from bobkit.buccaneer import Ca_sequence
  >>> Ca_sequence.set_cpus(1)
  >>> Ca_sequence.set_semet(False)
  >>> # sequence reliability used in Buccaneer is 0.95, default: 0.5
  >>> caseq = Ca_sequence(0.95)
  >>> caseq(mol_wrk, xwrk, llktargets, seq_wrk)
  >>> # return number of C-alphas sequenced
  >>> caseq.num_sequenced()
  >>> # return sequencing results up to top 5 scores for each chain as string
  >>> results = caseq.format()
  >>> print(results)


Sequence_score_threaded class
-----------------------------
.. autosummary::

   bobkit.buccaneer.Sequence_score_threaded

Sequence_threaded class
-----------------------
.. autosummary::

   bobkit.buccaneer.Sequence_threaded
