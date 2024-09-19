.. highlight: python

Growing chain fragments
=======================

The growing stage as it suggests will try to add C-alpha groups to either end of the C-alpha positions or input chains.
This is performed by optimising the log-likelihood fit to density for the new C-alpha groups while not disobeying 
to the constraints of the Ramachandran plot. The log-likelihood function used for evaluating C-alpha positions added
by growth is the same used in the initial finding stage. However, it is now evaluated in real space for each
candidate position and orientation instead of using the fast-fourier transform (FFT) approach.

Classes which are used in the growing routine:

* :ref:`Ca_grow <cagrow>` for growing C-alphas using density
* :ref:`Grow_threaded <growthread>` for growing C-alpha groups
* :ref:`Target_fn_refine_n_terminal_build <tgtref_nterm>` for refining grown C-alpha groups
* :ref:`Target_fn_refine_c_terminal_build <tgtref_cterm>` for refining grown C-alpha groups

.. _cagrow:

Ca_grow class
-------------

.. autosummary::

   bobkit.buccaneer.Ca_grow

The :py:meth:`constructor <bobkit.buccaneer.Ca_grow.__init__>` takes a maximum number 
of residues to add in either direction. The instance can then be :py:meth:`called <bobkit.buccaneer.Ca_grow.__call__>`
by giving it a :py:obj:`MiniMol <bobkit.clipper.MiniMol>` instance, a map and a target. 
It will first find starting chains to expand.
A cutoff threshold for the log-likelihood function is determined through map statistics. 
It is required to determine when to stop growing the chain in either direction. . 

Several optimisations are used to improve the performance. The log-likelihood function
is approximated by only using the grid points in the calculation.
The best 50 conformations of the first residue are used to build the second residue
and the best 30 combined scores are then rescored using all hte points in the
log-likehood function. The simplex algorithm search is used to refine the
Ramachandran angles for the best solution.

:py:meth:`Ca_grow.set_cpus() <bobkit.buccaneer.Ca_grow.set_cpus>` can be used to set 
the number of cpu threads to use in running the calculations.

.. doctest::

  >>> from bobkit.buccaneer import Ca_grow 
  >>> Ca_grow.set_cpus(1)
  >>> cagrow = Ca_grow(25)
  >>> # mol_wrk - MiniMol instance of working model 
  >>> # xwrk - working map
  >>> # llktgt - likelihood map target
  >>> cagrow(mol_wrk, xwrk, llktgt)

.. _growthread:

Grow_threaded class
-------------------
.. autosummary::

   bobkit.buccaneer.Grow_threaded
  
.. _tgtref_nterm:

Target_fn_refine_n_terminal_build class
---------------------------------------
.. autosummary::

   bobkit.buccaneer.Target_fn_refine_n_terminal_build

.. _tgtref_cterm:

Target_fn_refine_c_terminal_build class
---------------------------------------
.. autosummary::
   
   bobkit.buccaneer.Target_fn_refine_n_terminal_build