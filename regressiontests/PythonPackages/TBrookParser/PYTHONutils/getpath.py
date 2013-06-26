def getpath(primaryFile,secondaryFile):

    """
    Given the path of a primaryFile and the name of the secondaryFile
    ensure path of secondaryFile is correct

    Assumption: the secondaryFile lives in the same directory as the primaryFile
    """ 

    import string
    p = primaryFile.split('/')
    q = secondaryFile.split('/')
    L = p[:-1]
    L.append(q[-1])

    return string.join(L,'/')

if __name__ == '__main__':
    pfile = 'joe/jim/testprimary'
    sfile = 'bob/tom/testsecondary'
    B = getpath(pfile,sfile)
    print B

