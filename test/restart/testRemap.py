#!/usr/bin/env python3

import scipy as sp
import scipy.linalg as spla

import truchas

def test_remap(tenv):
    """Test creating and reading a restart file"""
    nfail = 0

    # run with 2, then 4 processors
    stdout, output2 = tenv.truchas(2, "ds11.inp", output_dir="remap_out2pe")
    stdout, output4 = tenv.truchas(4, "ds11.inp", output_dir="remap_out4pe")

    tol = 1e-2
    sid2 = output2.num_series()
    sid4 = output4.num_series()

    # compare the final temperature fields
    field2 = output2.field(sid2, "Z_TEMP")
    field4 = output4.field(sid4, "Z_TEMP")

    max = spla.norm(field4 - field2,  sp.inf)
    min = spla.norm(field4 - field2, -sp.inf)
    l2  = spla.norm(field4 - field2)

    status = "PASS" if max <= tol else "FAIL"
    print("{:s}: temperature: min {:8.2e} max {:8.2e} L2 {:8.2e}".format(status,min,max,l2))
    if max > tol: nfail += 1

    # compare the final temperature gradients
    field2 = output2.field(sid2, "Grad_T")
    field4 = output4.field(sid4, "Grad_T")

    max2 = spla.norm(field2,  sp.inf)
    max = spla.norm(field4 - field2,  sp.inf)
    min = spla.norm(field4 - field2, -sp.inf)
    l2  = spla.norm(field4 - field2)

    status = "PASS" if max/max2 <= tol else "FAIL"
    print("{:s}: grad_T: max of 2pe {:8.2e} min {:8.2e} max {:8.2e} L2 {:8.2e}".format(status,max2,min,max,l2))
    if max/max2 > tol: nfail += 1

    return nfail

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = test_remap(tenv)
    assert nfail == 0
