def getdir(filename):
    "Returns directory which filename lives in"
    import os
    return os.path.dirname(os.path.abspath(filename))
