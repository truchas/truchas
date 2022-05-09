.. _TOOLPATH_Namelist:

TOOLPATH Namelist
=================
The TOOLPATH namelist defines a path through space, especially that taken
by a machine tool, such as a laser or build platform, in the course of a
manufacturing process. Its use is not limited to such cases, however. The
path is specified using a simple command language that is adapted to common
CNC machine languages. The current implementation is limited to a path
through Cartesian 3-space (3-axis). The command language is described at
the end of the section.

.. note::

   :Required/Optional: Optional
   :Single/Multiple Instances: Multiple

Namelist Variables
--------------------------

.. contents::
   :local:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A unique name used to identify a particular instance of this namelist.
Clients will reference the toolpath using this name.

:Type: case sensitive string (31 characters max)
:Default: none


command_string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A string from which to read the commands that define the toolpath. This is
only suitable for very simple paths; use `command_file`_ instead for more
complex paths.

:Type: string (1000 characters max)
:Default: none

.. tip::
   Use single quotes to delimit the string to avoid conflicts with double
   quotes used within the commands.


command_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The path to a file from which to read the commands that define the toolpath.
If not an absolute path, it will be interpreted as a path relative to the
Truchas input file directory.

:Type: string
:Default: none

.. tip::
   C++ style comments may used in the file; all text from the "//"
   to the end of the line is ignored.


start_time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The starting time of the toolpath.

:Type: real
:Default: 0


start_coord
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The starting coordinates of the toolpath.

:Type: real 3-vector
:Default: (0,0,0)


time_scale_factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional multiplicative factor by which to scale all time values. This
is applied to all relevant namelist variables as well as toolpath commands,
including their speeds and accelerations.

:Type: real
:Default: 1
:Valid values: :math:`\gt 0`


coord_scale_factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional multiplicative factor by which to scale all time values. This
is applied to all relevant namelist variables as well as toolpath commands,
including their speeds and accelerations.

:Type: real
:Default: 1
:Valid values: :math:`> 0`


write_plotfile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Enable this flag to have a discrete version of the toolpath written to a
disk file. The file will be located in the Truchas output directory and
be named `toolpath-<name>.dat`, where `<name>` is the name assigned to the
namelist. If enabled, `plotfile_dt`_ must be specified.

:Type: logical
:Default: false

.. note::
    The file is a multi-column text file where each line gives the
    toolpath data at a specific time. The columns are, in order, the segment
    index, the time, the three position coordinates, and the flag settings
    (0 for clear and 1 for set). Not all flags are written, only those that
    were set at some point. The initial comment line starting with a "#"
    labels the columns. Data is written at the end point times of each path
    segment, and at zero or more equally spaced times within each segment
    interval. The latter frequency is determined by `plotfile_dt`_.


plotfile_dt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Output time frequency used when writing the toolpath to a disk file.
See `write_plotfile`_.

:Type: real
:Default: none
:Valid values: :math:`> 0`


partition_ds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Assign a value to this parameter to generate an additional discrete version
of the toolpath that is required by some clients. The value specifies the
desired spacing in path length. Refer to the client's documentation on
whether this is needed and for further details.

:Type: real
:Default: none
:Valid values: :math:`> 0`

.. note::
    Some clients require a discete version of the toolpath that consists
    of a time-ordered sequence of (time, coordinate) pairs. The sequence
    includes the end points of the segments. In addition each segment is
    partitioned into one or more parts of equal path length approximately
    equal to, but no greater than `partition_ds`_.


The Toolpath Command Language
---------------------------------
A toolpath is represented as a continuous sequence of :math:`2+n` path
segments defined on time intervals :math:`(-\infty,t_0]`, :math:`[t_0,t_1]`,
:math:`[t_1,t_2]`, ..., :math:`[t_n,\infty)`. The initial and final unbounded
path segments are automatically generated and are not specified. The
individual segments are simple paths (e.g., no motion or linear motion) that
are defined by path commands. The initial segment is a dwell at the position
given by `start_coord`_ that ends at time :math:`t_0` given by `start_time`_.
Likewise the final segment is a dwell starting at a time and at a position
determined by the preceding sequence of path segments. A set of flags
(0 through 31) is associated with each path segment. Clients of a toolpath
can use the setting of a flag (set or clear) for a variety of purposes; for
example, to indicate that a device is on or off for the duration of the
segment. All flags start clear in the initial segment. The toolpath command
language is expressed using JSON text. Except for the integer flag numbers,
the real numeric values may be entered as integer or floating point numbers;
that is, "1" is an acceptable alternative to "1.0". The specification of a
toolpath takes the form

    [ `command`, `command`, ... ]

where `command` is one of the following:

["dwell", :math:`dt` ]
    Remain at the current position for the time interval :math:`dt`.

["moverel", [ :math:`dx,dy,dz` ], :math:`s, a, d` ]
   Linear displacement from the current position. Motion accelerates from
   rest to a constant speed, and then decelerates to rest at the position
   :math:`(dx,dy,dz)` relative to the current position. The linear speed,
   acceleration, and deceleration are :math:`s`, :math:`a`, and :math:`d`,
   respectively. If :math:`d` is omitted, its value is taken to be :math:`a`.
   If both :math:`a` and :math:`d` are omitted then instantaneous acceleration
   (deceleration) to (from) speed :math:`s` is assumed.

["setflag", :math:`n_1, n_2,\ldots` ]
    Sets the listed flags. The setting of a flag holds for all subsequent
    motions (above) until changed.

["clrflag", :math:`n_1, n_2,\ldots` ]
    Clears the listed flags. The setting of a flag holds for all subsequent
    motions (above) until changed.
