.. highlight:: console

Running Buccaneer in Python
===========================

A python implementation of Buccaneer (buccaneer.py) can be found inside the **scripts** folder in the repository.

Built model requires further refinement as Buccaneer_ itself does not perform refinement of the model.

.. literalinclude:: buccaneer-help.txt
    :language: console

The usual routine is to run 2-3 cycles of Buccaneer_ and then refine the model before using it as a starting model for further building.
The build-refine cycle can be performed for several times in an attempt to complete the model.

The default library reference structure and data used by CCP4's Buccaneer_ hardly needs changing.
These reference data files can be found inside the **include/reference_data** folder in the repository.
Reference data 1TQW is used for X-ray cases, while EMD-4116 is used for EM cases.

.. _Buccaneer: https://www.ccp4.ac.uk/html/cbuccaneer.html

References
----------
1. `Cowtan, K. (2006) The Buccaneer software for automated model building. Acta Cryst. D62, 1002-1011`_

2. `Cowtan, K. (2010) Recent developments in classical density modification. Acta Cryst. D66, 470–478`_

3. `Cowtan, K. (2012) Completion of autobuilt protein models using a database of protein fragments. Acta Cryst. D68, 328-335`_

4. `Hoh, S.W., Burnley, T. & Cowtan, K. (2020) Current approaches for automated model building into cryo-EM maps using Buccaneer with CCP-EM. Acta Cryst. D76, 531-541`_

.. _Cowtan, K. (2006) The Buccaneer software for automated model building. Acta Cryst. D62, 1002-1011: https://doi.org/10.1107/S0907444906022116
.. _Cowtan, K. (2010) Recent developments in classical density modification. Acta Cryst. D66, 470–478: https://doi.org/10.1107/S090744490903947X
.. _Cowtan, K. (2012) Completion of autobuilt protein models using a database of protein fragments. Acta Cryst. D68, 328-335: https://doi.org/10.1107/S0907444911039655
.. _Hoh, S.W., Burnley, T. & Cowtan, K. (2020) Current approaches for automated model building into cryo-EM maps using Buccaneer with CCP-EM. Acta Cryst. D76, 531-541: https://doi.org/10.1107/S2059798320005513