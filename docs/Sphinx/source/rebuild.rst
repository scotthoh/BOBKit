.. highlight: python

Rebuilding side chains/atoms
============================

At this stage, side chain atoms and carbonyl oxygens are rebuilt.

Ca_build class
--------------
.. autosummary::

   bobkit.buccaneer.Ca_build

The :py:meth:`constructor <bobkit.buccaneer.Ca_build.__init__>` takes in two arguments:

* ``newrestype`` - three-letter code, new residue will be built as given type
* ``flexible`` - boolean type to set if rotamers are flexible

An instance of the class Ca_build can be :py:meth:`called <bobkit.buccaneer.Ca_build.__call__>`
by passing it the following arguments:

* ``mol`` - working model
* ``xmap`` - working map

.. doctest::

  >>> from bobkit.buccaneer import Ca_build
  >>> cabuild = Ca_build("ALA", False)
  >>> cabuild(mol_wrk, xwrk)

The equivalent static method :py:meth:`build <bobkit.buccaneer.Ca_build.build>` 
takes in four arguments:

* ``mol`` - working model
* ``xmap`` - working map
* ``newrestype`` - three-letter code, new residue will be built as given type (optional)
* ``flexible`` - boolean type to set if rotamers are flexible (optional)

.. doctest::

  >>> Ca_build.build(mol_wrk, xwrk, "ALA", False)
  >>> # or just the following where by default newrestype="ALA", flexible=False
  >>> Ca_build.build(mol_wrk, xwrk)