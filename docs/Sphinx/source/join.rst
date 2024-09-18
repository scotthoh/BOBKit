.. highlight: python

Joining chain fragments
=======================

The joining stage will try to merge the many overlapping chain fragments which may or may not be consistent with each other.

Classes which can be used for the growing routine:

* ``Ca_join`` for merging overlapped C-alpha chains
* ``Ca_join.Node`` to store path length and pointers to predecessor or successor.
* ``Ca_join.Tri_residue`` to store score of a fragment containing three residues

Ca_join class
-------------
.. autosummary::

   bobkit.buccaneer.Ca_join

The constructor takes in two arguments: merge and join radii. The default values are 2.0 Å for both.
Calling the Ca_join instance will run the joining function. Firstly, chain fragments are split into
a series of overlapping tri-residue fragments. Each tri-residue fragment overlaps its neighbour by two residues.
Multiple traces of the same chain fragment are then merged by combining any pair of tri-residues
for which all three C-alpha atoms match to within the given merging radius (default: 2.0 Å).
This is achieved by averaging all the coordinates of each main chain atom of each tri-residue. The result is a 
model in which multiple consistent traces of the same chain segment have been removed.

These tri-residues are examined for potential reassembling into chains. Any pair of tri-residue
for which the second and third C-alpha atoms of the first tri-residue matches, to within a given joining radius 
(default: 2.0 Å), the first and second C-alpha atoms of the second tri-residue will be considered as a potential join.
The correct routing is determined by taking the longest possible path or non-looped chain.

Ca_join.Node class
------------------
.. autosummary::

   bobkit.buccaneer.Ca_join.Node

.. doctest::

  >>> from bobkit.buccaneer import Ca_join
  >>> cajoin = Ca_join()
  >>> # pass in a working model yield from previous growing stage
  >>> cajoin(mol_wrk)

This class also has a static function which can be called to run the joining routine.
It takes the following arguments:

* ``mol`` - a MiniMol object containing the working model
* ``rmerge`` - merge radius
* ``rjoin`` - join radius
* ``com`` - ASU centre for output of results

.. doctest::

  >>> # com: ASU centre for output of results
  >>> # just an example here
  >>> com = clipper.Coord_orth(0.0,0.0,0.0)
  >>> Ca_join.join(mol_wrk, 2.0, 2.0, com)

Ca_join.Tri_residue class
-------------------------
.. autosummary::

   bobkit.buccaneer.Ca_join.Tri_residue

