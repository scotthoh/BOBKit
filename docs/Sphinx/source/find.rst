.. highlight: python

Finding C-alphas
================

This stage aims to locate a few very probably C-alpha positions in the electron-density map
using a search target based on a 4 Å sphere centred on the C-alpha atom. A six-dimensional
search in both position and orientation is performed and the highest scoring matches are 
assumed to be correct. These C-alpha positions will be used as seed points from which longer
chains will be grown. In subsequent cycles, the finding step is modified to preferentially
find C-alpha positions which are in regions where no model is present.

Classes which are used in the finding routine:

* :ref:`Ca_find <cafind>` for finding C-alphas in density
* :ref:`Search_threaded <searchthreaded>` for searching C-alpha groups
* :ref:`SearchResult <searchresult>` to store result scores.
* :ref:`SSfind <ssfind>` alternative to fffear for fast secondary structure finding
* :ref:`Target_fn_refine_llk_map_target <tgt_ref>` to refine C-alpha groups

.. _cafind:

Ca_find class
-------------
.. autosummary::

   bobkit.buccaneer.Ca_find

The :py:meth:`constructor <bobkit.buccaneer.Ca_find.__init__>` takes in two arguments.
An int for number of seeds to find and a float for resolution limit.

:py:meth:`Ca_find.set_cpus() <bobkit.buccaneer.Ca_find.set_cpus>` can be used to set the 
number of cpu threads to use in running the calculations.

There are two methods used in the finding routine. The fastest method typically gives 2-3x
speedup over the best method for a very similar model, but results might vary.

* :py:attr:`Ca_find.TYPE.SECSTRUC <bobkit.buccaneer.Ca_find.TYPE.SECSTRUC>` for fastest method
* :py:attr:`Ca_find.TYPE.LIKELIHOOD <bobkit.buccaneer.Ca_find.TYPE.LIKELIHOOD>` for best method

Either of the types can be passed as argument when :py:meth:`calling the Ca_find <bobkit.Ca_find.__call__>` 
instance. The following arguments are required to run the finding routine:-

* a **clipper.MiniMol** object containing the working model or can be empty (during the first find), constructed from the spacegroup and unit cell of the input data
* a **buccaneer.KnownStructure** object containing information of known structure (can be empty if no known structure)
* a working map of type **clipper.XMap_float**
* a likelihood map target of type **buccaneer.LLK_map_target**
* fastest or best finding method **Ca_find.TYPE.SECSTRUC** or **Ca_find.TYPE.LIKELIHOOD**
* modelindex: 0 to use all results, otherwise downweight 50% of results on the basis of position in ASU when filtering prior model (input work model)

The C-alpha finding stage starts by making a prior map from the given model and working map.
The prior map will be turned into a z-score followed by the filtering of the prior based on multi-model index given.
Search is performed and a long scores list is created from the maps of results. A pruned scores lists is then created omitting 
near-clashes. The best matches are refined before building a backbone model as the output.
The working model passed in to the finding method will be updated with C-alpha groups of the backbone model built.

.. doctest::

   >>> from bobkit.buccaneer import Ca_find
   >>> Ca_find.set_cpus(1)
   >>> cafind = Ca_find(500, 2.0)
   >>> cafind(mol_wrk, knownstruc, xwrk, llktgt, findtype, modelindex)

.. _searchthreaded:

Search_threaded class
---------------------
.. autosummary::

   bobkit.buccaneer.Search_threaded

.. _searchresult:

SearchResult class
------------------
.. autosummary::

   bobkit.buccaneer.SearchResult

.. _ssfind:

SSfind class
------------
.. autosummary::

   bobkit.buccaneer.SSfind

.. _tgt_ref:

Target_fn_refine_llk_map_target class
-------------------------------------
.. autosummary::

   bobkit.buccaneer.Target_fn_refine_llk_map_target