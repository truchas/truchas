.. _TOOLPATH_Namelist:

.. toctree::
   :maxdepth: 1

TOOLPATH Namelist (Experimental)
==================================

Overview
------------
The :ref:`TOOLPATH<TOOLPATH_Namelist>` namelist defines a path through space, especially that taken by a machine tool, such as a laser or build platform, in the course of a manufacturing process. Its use is not limited to such cases, however. The path is specified using a simple command language that is adapted to common CNC machine languages. The current implementation is limited to a path through Cartesian 3-space (3-axis). The command language is described at the end of the chapter.

TOOLPATH Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`Name<TP_Name>`
* :ref:`Command_String<TP_CS>`
* :ref:`Command_File<TP_CF>`
* :ref:`Start_Time<TP_ST>`
* :ref:`Start_Coord<TP_SC>`
* :ref:`Time_Scale_Factor<TP_TSF>`
* :ref:`Coord_Scale_Factor<TP_CSF>`
* :ref:`Write_Plotfile<TP_WP>`
* :ref:`Plotfile_Dt<TP_PD>`
* :ref:`Partition_Ds<TP_PDs>`

.. _TP_Name:

Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist. Clients will reference the toolpath using this name. 
| **Type**        : case sensitive string (31 characters max)
| **Default**     : none

.. _TP_CS:

Command_String
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A string from which to read the commands that define the toolpath. This is only suitable for relatively simple paths; use :ref:`Command_File<TP_CF>` instead for more complex paths.
| **Type**        : string (1000 characters max)
| **Default**     : none
| **Notes**       : Use single quotes to delimit the string to avoid conflicts with double quotes used within the commands.

.. _TP_CF:

Command_File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The path to a file from which to read the commands that define the toolpath. If not an absolute path, it will be interpreted as a path relative to the Truchas input file directory.
| **Type**        : string 
| **Default**     : none
| **Notes**       : C++ style comments may used in the file; all text from the ‘//’ to the end of the line is ignored.

.. _TP_ST:

Start_Time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The starting time of the toolpath.
| **Type**        : real 
| **Default**     : 0

.. _TP_SC:

Start_Coord
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The starting coordinates of the toolpath.
| **Type**        : real 3-vector
| **Default**     : (0,0,0)

.. _TP_TSF:

Time_Scale_Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An optional multiplicative factor by which to scale all time values. This applies to all namelist variables as well as toolpath commands.
| **Type**        : real 
| **Default**     : 1
| **Valid values**: :math:`\gt 0`
| **Notes**       : This is applied appropriately to speeds and accelerations in the toolpath commands. 

.. _TP_CSF:

Coord_Scale_Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An optional multiplicative factor by which to scale all time values. This applies to all namelist variables as well as toolpath commands.
| **Type**        : real 
| **Default**     : 1
| **Valid values**: > 0
| **Notes**       : This is applied appropriately to speeds and accelerations in the toolpath commands.

.. _TP_WP:

Write_Plotfile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Enable this flag to have a discrete version of the toolpath written to a disk file. The file will be located in the Truchas output directory and be named ``toolpath-name.dat``, where **name** is the name assigned to the namelist. If enabled, :ref:`Plotfile_Dt<TP_PD>` must be specified.
| **Type**        : logical 
| **Default**     : False
| **Notes**       : The file is a multi-column text file where each line gives the toolpath data at a specific time. The columns are, in order, the segment index, the time, the three position coordinates, and the flag settings (0 for clear and 1 for set). Not all flags are written, only those that were set at some point. The initial comment line starting with a ‘#’ labels the columns. Data is written at the end point times of each path segment, and at zero or more equally spaced times within each segment interval. The latter frequency is determined by :ref:`Plotfile_Dt<TP_PD>`.

.. _TP_PD:

Plotfile_Dt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Output time frequency used when writing the toolpath to a disk file. See :ref:`Write_Plotfile<TP_WP>`.
| **Type**        : real
| **Default**     : none
| **Valid values**: :math:`> 0`

.. _TP_PDs:

Partition_Ds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Assign a value to this parameter to generate an additional discrete version of the toolpath that is required by some clients. The value specifies the desired spacing in path length. Refer to the client’s documentation on whether this is needed and for further information.
| **Type**        : real
| **Default**     : none
| **Valid values**: :math:`> 0`
| **Note**        : Some clients require a discete version of the toolpath that consists of a time-ordered sequence of (time, coordinate) pairs. The sequence includes the end points of the segments. In addition each segment is partitioned into one or more parts of equal path length approximately equal to, but no greater than :ref:`Partition_Ds<TP_PDs>`.

The Toolpath Command language
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A toolpath is represented as a sequence of :math:`n + 2` continuous path segments defined on time intervals :math:`(-\infty,t_0], [t_0,t_1], [t_1,t_2], . . . ,[t_n,\infty).` The individual segments are simple paths (e.g., no motion or linear motion) that are defined by path commands. A set of flags (0 through 31) is associated with each path segment. Clients of a toolpath can use the setting of a flag (set or clear) for a variety of purposes; for example, to indicate that a device is on or off for the duration of the segment.
The toolpath command language is expressed using JSON text. The specification of a toolpath takes the form

:math:`[command, command, ......]`

where :math:`command` is one of the following:

**["dwell", dt]** - Remain at the current position for the time interval :math:`dt`.

**["moverel", [dx,dy,dz], s, a, d]** - Linear displacement from the current position. Motion accelerates from rest to a constant speed, and then decelerates to rest at the position (dx,dy,dz) relative to the current position. The linear speed, accleration, and deceleration are **s**,**a**, and **d**, respectively. If **d** is omitted, its value is taken to be **a**. If both **a** and **d** are omitted then instantaneous acceleration/deceleration to speed **s** is assumed.

**["setflag",n1,n2,...]**
**["clrflag",n1,n2,...]** - Sets or clears the listed flags. 32 flags (0 through 31) are available. The setting of a flag holds for all subsequent motions (above) until changed.

Except for the integer flag numbers, the real numeric values may be entered as integer or floating point numbers; that is, "1" is an acceptable alternative to "1.0". The initial and final unbounded path segments are automatically generated and are not specified. The initial segment is a dwell at the position given by :ref:`Start_Coord<TP_SC>` that ends at time :math:`t_0` given by :ref:`Start_Time<TP_ST>`. Likewise the final segment is a dwell starting at a time and at a position determined by the preceding sequence of path segments. All flags start clear in the initial segment.





