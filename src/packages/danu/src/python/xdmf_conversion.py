#===============================================================================
#
#  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
#  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
#  in the LICENSE file found in the top-level directory of this distribution.
#
#===============================================================================

import os
import sys
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

from numpy import array, size, where, unique, int32

import Danu


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' \
            + tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def dtype_extract(d):
    dstr = str(d)
    conversion_table = {
            "int8": ("Int", "1"),
            "int16": ("Int", "2"),
            "int32": ("Int", "4"),
            "int64": ("Int", "8"),
            "float32": ("Float", "4"),
            "float64": ("Float", "8"),
        }
    if dstr in conversion_table:
        return conversion_table[dstr]
    else:
        raise ValueError("dtype '%s' is not implemented yet" % dstr)

def main(filename, filename_out, filename_mesh):
    xdmf_top = Element("Xdmf", {
        "Version": "2.0",
        "xmlns:xi": "http://www.w3.org/2001/XInclude",
    })
    xdmf_domain = SubElement(xdmf_top, "Domain")
    xdmf_grid_parent = SubElement(xdmf_domain, "Grid", {
        "GridType": "Collection",
        "CollectionType": "Temporal",
    })

    f = Danu.Output(filename)
    sim = f.get_simulation("MAIN")

    # Read the mesh
    mesh = sim.open_mesh_link()
    nodes = array(mesh.coordinates()).T
    elements = mesh.read_connectivity()
    #elements_offset = elements.attrs["Offset"][0]
    elements_offset = 1 # FIXME: Obtain this from Danu
    # FIXME: The path should be: "/Simulations/MAIN/Mesh/Element Connectivity", but
    # nested group creation is not implemented in the C code yet, so we use _
    # instead of / for now. Once this is fixed, also update the path below in
    # XML.
    h = Danu.danu_h5_create_handle()
    Danu.danu_h5_open(h, filename_mesh, "/Simulations_MAIN_Mesh")
    elements1d_size = Danu.danu_hex_save(h, elements, "Element Connectivity")

    try:
        blockid = sim.data_read('BLOCKID')
    except RuntimeError:
        print 'BLOCKID not found'
        raise
    ids = unique(blockid)
    element_blocks = [array(where(blockid == id)[0], dtype=int32) for id in ids]
    for i in range(len(ids)):
        Danu.danu_h5_save_int_1d_array(h, "Element Blocks %d" % i,
                element_blocks[i])
    Danu.danu_h5_close(h)
    Danu.danu_h5_free_handle(h)

    # Series
    for n in range(1, sim.sequence_count()+1):
        xdmf_grid = SubElement(xdmf_grid_parent, "Grid", {
            "GridType": "Uniform",
        })

        # Nodes
        xdmf_geometry = SubElement(xdmf_grid, "Geometry", {"GeometryType": "XYZ"})
        data_type, precision = dtype_extract(nodes.dtype)
        xdmf_dataitem = SubElement(xdmf_geometry, "DataItem", {
            "DataType": data_type,
            "Precision": precision,
            "Dimensions": "%d %d" % nodes.shape,
            "Format": "HDF",
            "Name": "Coordinates",
        })
        xdmf_dataitem.text="%s:/Simulations/MAIN/Mesh/Nodal Coordinates" % filename

        # Elements
        xdmf_topology = SubElement(xdmf_grid, "Topology", {
            "BaseOffset": str(elements_offset),
            "NumberOfElements": str(elements.shape[0]),
            "TopologyType": "Mixed",
        })
        data_type, precision = dtype_extract(elements.dtype)
        xdmf_dataitem = SubElement(xdmf_topology, "DataItem", {
            "DataType": data_type,
            "Precision": precision,
            "Dimensions": "%d" % elements1d_size,
            "Format": "HDF",
            "Name": "Connectivity",
        })
        xdmf_dataitem.text="%s:/Simulations_MAIN_Mesh/Element Connectivity" % filename_mesh

        # Sets
        for i in range(len(ids)):
            xdmf_set = SubElement(xdmf_grid, "Set", {
                "SetType": "Cell",
                "Name": "Block %d" % ids[i],
            })
            xdmf_dataitem = SubElement(xdmf_set, "DataItem", {
                "NumberType": "Int",
                "Dimensions": "%d" % size(element_blocks[i]),
                "Format": "HDF",
            })
            xdmf_dataitem.text="%s:/Simulations_MAIN_Mesh/Element Blocks %d" \
                    % (filename_mesh, i)


        # Time step
        seq_name = sim.get_sequence_name(n)
        time_step = sim.get_sequence(seq_name)
        xdmf_time = SubElement(xdmf_grid, "Time", {
            "TimeType": "Single",
            "Value": str(time_step.time),
        })

        # Loop over all fields in a given time step
        for field_name in time_step.data_list():
            try:
                attrs = time_step.data_attributes(field_name)
            except:
                attrs = {}
            dims = time_step.get_data_dimensions(field_name)
            if len(dims) == 1:
                type_str = "Scalar"
                dimensions_str = "%d" % dims[0]
            elif len(dims) == 2:
                if dims[1] == 3 and field_name in ["Z_VC", "Displacement"]:
                    # ParaView shows "Vector" type as having components (X, Y,
                    # Z) and can also show magnitude. We only want to use this
                    # for velocity.
                    type_str = "Vector"
                else:
                    type_str = "Matrix"
                dimensions_str = "%d %d" % tuple(dims)
            else:
                raise ValueError("Fields with more than 2 dimensions are not supported")
            # If the HDF5 field name is in this dictionary, the new name will be
            # used instead (e.g. in Paraview). Otherwise the name will remain
            # unchanged.
            field_names = {
                    "Z_VC": "Velocity",
                }

            if dims[0] == nodes.shape[0]:
                cell_type = "Node"
            elif dims[0] == elements.shape[0]:
                cell_type = "Cell"
            else:
                raise Exception("Unsupported array size (neither Node nor Cell).")

            if type_str == "Matrix":
                # XDMF type 'Matrix' should work for a matrix (p, m), but for
                # some reason Paraview cannot show the m components, so we
                # instead save 'm' scalar arrays of dimension 'p'.
                for m in range(dims[1]):
                    output_field_name = attrs.get("FIELDNAME%d" % (m+1))
                    if not output_field_name:
                        output_field_name = field_names.get(field_name,
                                field_name) + " (Field %d)" % (m+1)
                    xdmf_attribute = SubElement(xdmf_grid, "Attribute", {
                        "Center": cell_type,
                        "Name": output_field_name,
                        "Type": "Scalar",
                    })
                    hyper_slab = SubElement(xdmf_attribute, "DataItem", {
                        "ItemType": "HyperSlab",
                        "Dimensions": "%d 1" % dims[0],
                        "Type": "HyperSlab",
                        "Name": "Slab",
                    })

                    xdmf_dataitem = SubElement(hyper_slab, "DataItem", {
                        "Dimensions": "3 2",
                        "Format": "XML",
                    })
                    xdmf_dataitem.text="0 %d\n1 1\n%d 1" % (m, dims[0])

                    field_dtype = "float64" # FIXME: Obtain this from Danu
                    data_type, precision = dtype_extract(field_dtype)
                    xdmf_dataitem = SubElement(hyper_slab, "DataItem", {
                        "DataType": data_type,
                        "Precision": precision,
                        "Dimensions": dimensions_str,
                        "Format": "HDF",
                        "Name": field_name,
                    })
                    xdmf_dataitem.text="%s:/Simulations/MAIN/Series Data/Series %d/%s" \
                            % (filename, n, field_name)
            else:
                output_field_name = attrs.get("FIELDNAME")
                if not output_field_name:
                    output_field_name = field_names.get(field_name, field_name)
                xdmf_attribute = SubElement(xdmf_grid, "Attribute", {
                    "Center": cell_type,
                    "Name": output_field_name,
                    "Type": type_str,
                })
                field_dtype = "float64" # FIXME: Obtain this from Danu
                data_type, precision = dtype_extract(field_dtype)
                xdmf_dataitem = SubElement(xdmf_attribute, "DataItem", {
                    "DataType": data_type,
                    "Precision": precision,
                    "Dimensions": dimensions_str,
                    "Format": "HDF",
                    "Name": field_name,
                })
                xdmf_dataitem.text="%s:/Simulations/MAIN/Series Data/Series %d/%s" \
                        % (filename, n, field_name)

    open(filename_out, "w").write(prettify(xdmf_top))

def start():
    if len(sys.argv) != 2:
        print 'Usage: xdmf-parser.py DANU_FILE.h5'
        sys.exit(1)
    filename = sys.argv[1]
    filename_out = os.path.splitext(filename)[0] + ".xmf"
    filename_mesh = os.path.splitext(filename)[0] + "_mesh" + ".h5"
    main(filename, filename_out, filename_mesh)
