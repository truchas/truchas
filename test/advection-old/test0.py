#!/usr/bin/env python3

import truchas

# Really nothing we can test here, other than it ran successfully
if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    stdout, output = tenv.truchas(4, "advection-old-0.inp")
