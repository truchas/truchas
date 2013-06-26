def unique(l):
    "Return a list of the elements in s, but without duplicates."

    from sets import Set
    u = []
    try:
        l2 = l.tolist()
    except:
        l2 = l
    u.extend(Set(l2))
    return u

if __name__ == '__main__':
    L = ['a', 'b', 'c', 'a', 'b', 'c']
    print L
    B = unique(L)
    print B
