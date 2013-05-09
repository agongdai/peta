def cmp(v1, v2):
    if v1 == v2:
        return 0
    if v1 > v2:
        return 1
    return -1

def slist_binary(sl, q):
    start = 0
    end = len(sl) - 1
    middle = 0
    cmp_rs = 0
    while (start <= end):
        middle = (start + end) / 2
        cmp_rs = cmp(q, sl[middle])
        if (cmp_rs == 0):
            return middle
        if (cmp_rs > 0):
            start = middle + 1
        else:
            end = middle - 1
    return -1

def slist_ins_pos(sl, q):
    start = 0
    end = len(sl) - 1
    middle = 0
    cmp_rs = 0
    while (start <= end):
        middle = (start + end) / 2
        cmp_rs = cmp(q, sl[middle])
        if (cmp_rs == 0):
            return middle
        if (cmp_rs > 0):
            start = middle + 1
        else:
            end = middle - 1
    if (end < 0):
        return -1
    if (start > len(sl) -1):
        return len(sl)
    return start

def rev_comp(sequence):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

print rev_comp('GTATCATTGAGCCCTCTCTAAAGGC')

#sl = [0, 1, 2.5, 3, 4, 5, 6, 7]
#print 0, slist_binary(sl, 0)
#print 1, slist_binary(sl, 1)
#print 2, slist_binary(sl, 2)
#print 3, slist_binary(sl, 3)
#print 4, slist_binary(sl, 4)
#print 5, slist_binary(sl, 5)
#print 6, slist_binary(sl, 6)
#print 7, slist_binary(sl, 7)
#
#print 8, slist_binary(sl, 8)
#print -1, slist_binary(sl, -1)
#print 0.5, slist_binary(sl, 0.5)
#print 1.5, slist_binary(sl, 1.5)
#print 2.5, slist_binary(sl, 2.5)
#print 3.5, slist_binary(sl, 3.5)
#print 4.5, slist_binary(sl, 4.5)
#print 5.5, slist_binary(sl, 5.5)
#print 6.5, slist_binary(sl, 6.5)
#print 7.5, slist_binary(sl, 7.5)
#print 8.5, slist_binary(sl, 8.5)