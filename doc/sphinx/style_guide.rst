.. role:: rst(code)
   :language: RST

Sphinx Documentation Style Guide
================================

This style guide provides best-practices for writing the Truchas Sphinx documentation. We encourage
documentation authors to follow this guide to ensure the quality and consistency of all Truchas
documentation.

.. contents:: Contents
   :local:
   :backlinks: none



Authors
-------
Every RST document should include the author near the top of the document. This is done with the
`sectionauthor directive
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-sectionauthor>`_.

.. code-block:: RST

   .. sectionauthor:: name <email>



Headings
--------
reStructuredText does not assign heading levels to particular characters. Instead, heading levels
are determined by the order in which they appear. For example, the first heading style creates a top
level heading, and any further instances of that heading style also creates top level heading.

For consistency, we use the following succession of headings.

.. code-block:: RST

   Page Title
   ==========

   Heading 1
   ---------

   Heading 2
   +++++++++

   Heading 3
   ^^^^^^^^^

   Heading 4
   ~~~~~~~~~



Figures
-------

Single Figure
+++++++++++++
When placing a single centered figure with no text wrapping around it, use a :rst:`figwidth` of at
most 90%. This provides a visual distinction between the text and the figure. The following code
block produces the figure below it.

.. code-block:: RST

   .. figure:: https://via.placeholder.com/800x600
      :figwidth: 90%
      :align: center

      This is a caption for the figure.

.. figure:: https://via.placeholder.com/800x600
   :figwidth: 90%
   :align: center

   This is a caption for the figure.

Notice how the figure has some padding to its right and left that distinguishes it from the text
before and after it. This padding is not part of the CSS to prevent excessive white space
around figures that have text wrapping around them.

Side-by-Side Figures
++++++++++++++++++++
Placing figures side-by-side in reStructuredText can be tricky. The best method we've found so far
is to put the figures in a table. We use the :rst:`list-table` directive when each figure has its
own caption, and the generic :rst:`table` directive when a caption is shared between figures.

Individual Captions
^^^^^^^^^^^^^^^^^^^
The :rst:`list-table` directive is useful when each figure has its own caption. For example, the
following code produces the figures below it.

.. code-block:: RST

   .. list-table::
      :align: center

      * - .. figure:: https://via.placeholder.com/400x300
             :width: 100%
             :align: center

             This is the caption.

        - .. figure:: https://via.placeholder.com/400x300
             :width: 100%
             :align: center

             This is the caption.

.. list-table::
   :align: center

   * - .. figure:: https://via.placeholder.com/400x300
          :width: 100%
          :align: center

          This is the caption.

     - .. figure:: https://via.placeholder.com/400x300
          :width: 100%
          :align: center

          This is the caption.

Shared Captions
^^^^^^^^^^^^^^^
The :rst:`list-table` syntax does not allow multi-column cells, so we instead use the more general
:rst:`table` directive. Due to how tables are parsed, however, this means that images must be
included indirectly through a substitution.

The following code block shows how to implement a caption spanning multiple images.

.. code-block:: RST

   .. |image1| image:: https://via.placeholder.com/300x500
      :width: 100%
      :align: middle

   .. |image2| image:: https://via.placeholder.com/300x500
      :width: 100%
      :align: middle

   .. |image3| image:: https://via.placeholder.com/300x500
      :width: 100%
      :align: middle

   .. table::
      :align: center

      +------------+------------+------------+
      |  |image1|  |  |image2|  |  |image3|  |
      +------------+------------+------------+
      | This caption spans all the images.   |
      | It describes the first image (left), |
      | the second image (center), and the   |
      | third image (right).                 |
      +--------------------------------------+

.. |image1| image:: https://via.placeholder.com/300x500
   :width: 100%
   :align: middle

.. |image2| image:: https://via.placeholder.com/300x500
   :width: 100%
   :align: middle

.. |image3| image:: https://via.placeholder.com/300x500
   :width: 100%
   :align: middle

.. table::
   :align: center

   +------------+------------+------------+
   |  |image1|  |  |image2|  |  |image3|  |
   +------------+------------+------------+
   | This caption spans all the images.   |
   | It describes the first image (left), |
   | the second image (center), and the   |
   | third image (right).                 |
   +--------------------------------------+



Namelists
---------
When describing a namelist and its parameters, use the following format. The exact heading level is
not important, as long as the parameters are sub-headings of the namelist.

.. code-block:: RST

   EXAMPLE Namelist
   ----------------

   PARAMETER_1
   +++++++++++
   Short description of parameter 1.

   .. namelist_parameter::
      :type: INTEGER
      :domain: parameter_1 >= 0
      :default: parameter_1 = 1

   A more in-depth description of parameter 1. This can include any content
   (e.g. images, math, code blocks) and is left to the author's discretion.

For examples of how namelists are described using this format, see the :doc:`PATCHES namelist
documentation <tools/RadE/patches/patches_namelist>`.

Namelist Parameter Directive
++++++++++++++++++++++++++++
When describing a namelist parameter, the 'namelist_parameter' directive should be used after the
short description. This custom directive is provided by the Truchas Sphinx build.

The directive inserts an element of the *nml-param* class, ensuring that all namelist
parameters are styled consistently throughout the documentation.

The directive takes three mandatory options:

#. **:type:** the Fortran type of the namelist parameter
#. **:domain:** the set of valid values for the parameter
#. **:default:** the default value if the parameter is not specified

The following code block produces the namelist parameter element below it.

.. code-block:: RST

   .. namelist_parameter::
      :type: INTEGER
      :domain: parameter_1 >= 0
      :default: parameter_1 = 1

.. namelist_parameter::
   :type: INTEGER
   :domain: parameter_1 >= 0
   :default: parameter_1 = 1



Inline Fortran Code
-------------------
The Truchas Sphinx build provides a special 'fortran' role for writing inline Fortran code. For
example, the following RST code

.. code-block:: RST

   :fortran:`integer :: x`

produces the inline Fortran code :fortran:`integer :: x`.



Other Tips
----------

Math default role
+++++++++++++++++
For documents with a lot of inline math, you can set the default role to math by adding the
line

.. code-block:: RST

   .. default-role:: math

to the top of the document. Inline math can then be delimited by backticks (`) without explicitly
specifying the *:math:* role. For example, the line

.. code-block:: RST

   The sum of the first :math:`n` natural numbers is given by :math:`n(n+1)/2`.

simplifies to

.. code-block:: RST

   The sum of the first `n` natural numbers is given by `n(n+1)/2`.
