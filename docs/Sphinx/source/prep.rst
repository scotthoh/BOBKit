.. highlight: python

Preparing log-likelihood map matching targets
=============================================

The log-likelihood (LLK) map matching targets are prepared prior to running the ten stages of Buccaneer.
They are used mainly in finding C-alphas positions, growing fragments, linking fragments, sequencing, 
and correcting insertions and deletions.

LLK_map_target class
--------------------
.. autosummary::
  
   bobkit.buccaneer.LLK_map_target

In Buccaneer, this class is used in determining the LLK fit of some desired density features from some 
region of a noisy electron density map. It contains:-

* methods to accumulate the LLK target from a number of sample density
  regions from a map with similar noise levels to the target map
* methods for a FFFear 6-dimensional (rotation/orientation) search of the target map
* methods for testing indivdual sample positions and orientations.

.. note::
  The results from the 6-d search and the fast and full LLK
  calculations are on different scales and so cannot be compared

The constructor takes in three arguments:

* ``rad`` - radius, in Angstroms (Å), in which LLK function will be used
* ``sampling`` - sampling spacing, in Angstroms (Å), for the LLK function
* ``type`` - sampling type: normal or correlation (optional, default: NORMAL)
  
  * ``LLK_map_target::TYPE::NORMAL`` or ``LLK_map_target::TYPE::CORREL``

``LLK_map_target.accummulate()`` is used to accumulate density statistics from 
a sample orientation in a sample map. This method takes in two arguments:

* ``xmap`` - sample map
* ``rtop`` - rotation-translation operator that generates the target in question

``LLK_map_target.prep_llk()`` is used to prepare LLK target after accumulating
density or loading map.

.. doctest::

  >>> from bobkit.buccaneer import LLK_map_target, Ca_group
  >>> # example of trying to get main chain targets
  >>> # mmol is of class clipper.MiniMol
  >>> # get residue of any type besides glycine and proline
  >>> res1 = mmol[0][1]
  >>> ca = Ca_group(res1)
  >>> rtop = ca.rtop_from_std_ori()
  >>> llktgt = LLK_map_target(4.0, 0.5, LLK_map_target.NORMAL)
  >>> llktgt.accumulate(xmap, rtop)
  >>> # prepare LLK target after accumulating density
  >>> llktgt.prep_llk()
  >>> # generate distribution of LLK values for a given map
  >>> llktgt.prep_llk_distribution(xmap_wrk)


The main chain and side chain LLK map targets can be prepared with the help of the
:ref:`ca_prep_ref` class. 


LLK_TargetList class
--------------------
.. autosummary:: 

   bobkit.buccaneer.LLK_TargetList

This class is used to hold a list of LLK_map_target. In buccaneer calculation, 
it is mainly used to hold the main chain and twenty side chain map targets.
The constructor can take in an existing list of LLK_map_target or an integer of the size
of vector to initialise.

.. doctest::

  >>> from bobkit.buccaneer import LLK_TargetList, LLK_map_target
  >>> # construct from existing list of LLK_map_target
  >>> list_llk = [LLK_map_target] * 20
  >>> llkcls = LLK_TargetList(list_llk)
  >>> # construct from size for list
  >>> llkcls = LLK_TargetList(20)
  >>> # to return a reference to the list
  >>> llkcls.get_vector()


.. _ca_prep_ref:

Ca_prep class
-------------
.. autosummary::

   bobkit.buccaneer.Ca_prep

This class is used for preparing main chain (C-alphas) and side chain (C-beta) 
log-likelihood map targets. In Buccaneer calculation, a radius of 4.0 Å is used 
for the main chain LLK function, while the side chain LLK function uses a radius 
of 5.5 Å.

The constructor takes in six arguments:

* ``main_tgt_rad`` - main chain target radius, in Angstroms (Å), to be used in LLK function
* ``side_tgt_rad`` - side chain target radius, in Angstroms (Å), to be used in LLK function
* ``rama_flt`` - ramachandran filter data: all, helix, strand, or nonhelix
  
  * ``Ca_prep.Rama_flt.rama_flt_all``, ``Ca_prep.Rama_flt.rama_flt_helix``, ``Ca_prep.Rama_flt.rama_flt_strand``, ``Ca_prep.Rama_flt.rama_flt_nonhelix``

* ``correl`` - flag to use correlation for LLK map targets
* ``seqnc`` - flag to indicate whether sequencing step is turned on
* ``debug`` - flag to turn on/off debugging mode.

An instance to the class can be called by passing it the following:

* ``llktgt`` - LLK map target to hold the data for C-alphas
* ``llkcls`` - a list of LLK map targets for side chain targets
* ``mol`` - reference model
* ``xmap`` - reference density map of the reference model

The example steps to prepare the map targets and then use them for
sequencing are as follow:

.. doctest::

    >>> from bobkit.buccaneer import Ca_prep, Ca_sequence
    >>> from bobkit.buccaneer import LLK_map_target, LLK_TargetList
    >>> 
    >>> llktgt = LLK_map_target()    # for C-alphas
    >>> llkcls = LLK_TargetList(20)  # for all 20 amino acids
    >>> rama_fltr = Ca_prep.Rama_flt.rama_flt_all  # default
    >>> buccaneer.Ca_prep.set_cpus(1)
    >>> caprep = buccaneer.Ca_prep(4.0, 5.5, rama_fltr, True, True, False)
    >>> # mol_ref and xref are reference model and density map
    >>> caprep(llktgt, llkcls, mol_ref, xref)
    >>> # example of using the LLK map targets for sequencing
    >>> Ca_sequence.set_cpus(1)
    >>> Ca_sequence.set_semet(False)
    >>> caseq = Ca_sequence(0.95)
    >>> caseq(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)