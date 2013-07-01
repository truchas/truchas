import sys

def SystemExit(fpwatch=sys.stdout):

    "Procedure to exit completely and cleanly from the postprocessor"

    print >> fpwatch
    print >> fpwatch, "Will not continue to postprocess this particular file."
    print >> fpwatch
    sys.exit()

