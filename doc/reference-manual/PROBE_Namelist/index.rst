.. _PROBE_Namelist:

.. toctree::
   :maxdepth: 1

PROBE Namelist
=============================

Overview
------------
The **PROBE** namelist is used to define the location in the computational domain where the value of specific solution quantities will be recorded at every time step. The data is written to a standard multi-column text file specified by the namelist. The solution time is written to the first column and the specified solution quantities to the remaining columns. Useful metadata about the probe is written to the first fewlines of the file. These lines begin with a # character and would be treated as comment lines by many post-processors. As many probes as desired may be defined, but each must write to a different file.

PROBE Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`coord<Prb_c>`
* :ref:`coord_scale_factor<Prb_csf>`
* :ref:`data<Prb_d>`
* :ref:`data_file<Prb_df>`
* :ref:`description<Prb_desc>`
* :ref:`digits<Prb_dgts>`


.. _Prb_c:

coord
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The spatial coordinates of the location of the probe. These coordinates may be further scaled by an optional scaling factor; see :ref:`coord_scale_factor<Prb_csf>` 
| **Type**        : real 3-vector
| **Default**     : none
| **Notes**       : For cell-centered quantities, data will be taken from the cell whose centroid is nearest this location, and for node-centered quantities data will be taken from the nearest node. 

.. _Prb_csf:

coord_scale_factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A multiplicative scaling factor applied to the coordinates of the probe.
| **Type**        : real 
| **Default**     : 1.0

.. _Prb_d:

data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The data quantity whose value will be recorded. The available options are:

   * "temperature - Cell-centered temperature
   * "pressure"   - Cell-centered fluid pressure
   * "velocity"   - Cell-centered fluid velocity
   * "vol-frac-*name* - Cell-centered volume fraction of material or phase *name*. Name may also be the reserved name "VOID".

| **Type**        : string 
| **Default**     : none

.. _Prb_df:

data_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of the probe output file. The file will be created in the output directory, and any existing file will be overwritten. Each probe must write to a different output file.
| **Type**        : string 
| **Default**     : none

.. _Prb_desc:

description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An optional text string that will be written to the header of the output file.
| **Type**        : string 
| **Default**     : none

.. _Prb_dgts:

digits
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The number of significant digits in the output data.
| **Type**        : integer 
| **Default**     : 6
