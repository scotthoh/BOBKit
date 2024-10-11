.. highlight: python

Linking chain fragments
=======================

At this stage, nearby N and C termini of every fragments are examined to see if they can be linked by a short loop.
The linking routine is done without regard to sequence. Thus, sometimes these linkages are made incorrectly but is
not serious as these resulting chains can be docked to the sequence. If the linkages are erroneous, the two parts of
the chain will usually dock to different places in the seuqence, which then can be corrected by breaking the chain again.

However, caution is needed when linking chains with previously assigned sequence. In the Buccaneer calculation, the use 
of the sequence to validate the link is no longer available when the chains are linked with previous assigned sequence.
So, it is necessary to limit the length of the bridging fragment to six amino acids (i.e. two amino acids overlap with 
each chain and a maximum of two amino acids of gap). Once the model is approaching completion, at which point errors 
become less likely (wrongly sequence fragments are morely likely to occur when fragments are short), the constraint to 
a maximum of six amino acids bridging fragment should probably be relaxed to allow longer loops to be built.

Classes which can be used for the linking routine:

* :ref:`Ca_link <calink>` for linking of nearby N and C termini of fragments
* :ref:`Ca_group <cagroup>` class to store and generate groups of N, C-alpha, and C atoms of neighour residues
* :ref:`ProteinLoop <proteinloop>` builder class for loops of various lengths

.. _calink:

Ca_link class
-------------
.. autosummary::

   bobkit.buccaneer.Ca_link
  
The :py:meth:`constructor <bobkit.buccaneer.Ca_link.__init__>` takes in two arguments: link radius (default: 5.0 Å, Buccaneer uses 10.0 Å) 
and torsion sampling (default: 24). An instance of the class can be :py:meth:`called <bobkit.buccaneer.Ca_link.__call__>`
by passing it the model/fragments to link, a working map and likelihood map target.

.. doctest::

  >>> from bobkit.buccaneer import Ca_link
  >>> calink = Ca_link(10.0, 24)
  >>> calink(mol_wrk, xwrk, llk_map_target)
  >>> # to check numbers of C-alphas linked
  >>> calink.num_linked()

.. _cagroup:

Ca_group class
--------------
.. autosummary::

   bobkit.buccaneer.Ca_group


.. _proteinloop:

ProteinLoop class
-----------------
.. autosummary::

   bobkit.buccaneer.ProteinLoop