Introduction
============

This is the introduction to the Glacier_Centerlines project.

**Motivation**
This project was created as a consequence of the increasing complexity of OGGM. Some people wanted to use a part of the features of the model without having to run it all. The feature that is wanted to be computable independentily from OGGM is glacier centerlines. Consequentily here we fullow the same methodology as in OGGM and reuse part of the code from OGGM.

.. _OGGM: https://OGGM.org
.. _OGGMgit: https://github.com/OGGM

There are 4 main scripts in the repository, listed below:

**Code organization**

.. admonition:: **Scripts**
    :class: info

    main.py
      The workflow of the project is executed in this script. It calls all the other scripts, loads DEM and glaciers outlines and computes the glacier heads, tails and ultimately the centerlines. 

    utils.py
      It is a Toolbox for the whole project. It provides some small functions and classes.

    functions.py
      Functions used in main.py. Most of the functions come from OGGM/core and OGGM/utils. Some others are new.
      
    params.py
      Here the main parameters used in main.py are listed and initialized.


