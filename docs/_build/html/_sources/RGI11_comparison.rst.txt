.. _RGI11_comparison:

RGI11-check
=======================================

(Corresponding to dev-RGI11 branch : https://github.com/GLIMS-RGI/glacier_centerlines/tree/dev-RGI11)

Results:
~~~~~~~~~~~
**Highlights:**

 - 84% out of 3927 were able to be computed (627 failed).

 - 2 results have been issued: main flowline and "all" flowlines (main + tributaries)
 	* Main flowline: most of the time new flowline matches OGGM. "Matches" means that it starts and ends at the same point, but with slightly different trajectories.
 	* All flowlines: Agreement is worse than in the main flowline results. Normally we get some more tributaries in our new tool.

**How has it been run?**

* Data
	- DEM taken from each glacier from RGI6-11. See a sample in `/test_data <https://github.com/GLIMS-RGI/glacier_centerlines/blob/dev-RGI11/test_data/>`_
	- Glacier outlines from RGI6-11 outlines. 
* Main script
	- `RGI11check.py <https://github.com/GLIMS-RGI/glacier_centerlines/blob/dev-RGI11/RGI11check.py>`_

**Detailed results:**

In purple: new results.

Dashed: OGGM old flowlines.
	
.. figure:: _static/RGI60-11.02209.png
   :width: 60%
   :align: left
	
   EX1 main flowline, differences determining main flowline: RGI60-11.02209


.. figure:: _static/RGI60-11.02337.png
   :width: 60%
   :align: left
   
   EX2 main flowline, differences determining heads and tails: RGI60-11.02337	



.. figure:: _static/RGI60-11.02242.png
   :width: 60%
   :align: left
   
   EX3 multiple flowlines, differences in amount of tributaries: RGI60-11.02242
   
.. figure:: _static/RGI60-11.02245.png
   :width: 60%
   :align: left   
   
   EX3 multiple flowlines, differences in amount of tributaries: RGI60-11.02245
   


Further info
------------

Check for more results yourselves in the data folder, main_flowlines.zip and multiple_flowlines.zip : https://github.com/GLIMS-RGI/glacier_centerlines/tree/dev-RGI11/test_data/
