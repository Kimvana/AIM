# local lib imports
from . import DataConversion as AIM_DC


# If any function here takes up too much time.... Is it worth it to switch to
# the python sets? These are, however, unordered....

def union_KvA(atlist1, atlist2):
    """
    Takes two sorted lists, and returns a sorted lists containing all elements
    that were present in any of the input lists.
    """
    Nat1 = len(atlist1)
    Nat2 = len(atlist2)

    at1c = 0
    at2c = 0

    unatlist = []

    if Nat1 == 0:
        if Nat2 == 0:
            pass
        else:
            unatlist = atlist2
    else:
        if Nat2 == 0:
            unatlist = atlist1
        else:
            while True:
                if atlist1[at1c] > atlist2[at2c]:
                    unatlist.append(atlist2[at2c])
                    at2c += 1
                elif atlist1[at1c] < atlist2[at2c]:
                    unatlist.append(atlist1[at1c])
                    at1c += 1
                else:
                    unatlist.append(atlist1[at1c])
                    at1c += 1
                    at2c += 1

                if at1c == Nat1 and at2c == Nat2:
                    break
                elif at1c == Nat1:
                    for i in range(at2c, Nat2, 1):
                        unatlist.append(atlist2[i])
                    break
                elif at2c == Nat2:
                    for i in range(at1c, Nat1, 1):
                        unatlist.append(atlist1[i])
                    break
            unatlist = unatlist

    return unatlist


def difference_KvA(atlist1, atlist2):
    """
    Takes two sorted lists and returns a sorted list containing all the
    elements that are present in atlist1, but not in atlist2
    """
    Nat1 = len(atlist1)
    Nat2 = len(atlist2)

    at1c = 0
    at2c = 0

    while True:
        if atlist1[at1c] > atlist2[at2c]:
            at2c += 1
        elif atlist1[at1c] < atlist2[at2c]:
            at1c += 1
        else:
            del atlist1[at1c]
            at2c += 1

        if at1c == Nat1:
            break
        elif at2c == Nat2:
            break

    return atlist1


def difference_c(atlist1, atlist2, clib):
    """
    The c way of calculating the difference of two lists (like difference_KvA)
    """
    lengths = [len(atlist1), len(atlist2), len(atlist1)]
    lengths_c = AIM_DC.ctypelist(lengths, 'int32')
    clib.group_difference(atlist1, atlist2, lengths_c)
    return atlist1


def threesort_KvA(atlist):
    """
    Sorts a list of three items.
    """
    if atlist[0] < atlist[1]:
        if atlist[0] < atlist[2]:
            out = [atlist[0]]
            if atlist[1] < atlist[2]:
                out.extend(atlist[1:])
            else:
                out.extend((atlist[2], atlist[1]))
        else:
            out = [atlist[2], atlist[0], atlist[1]]
    else:
        if atlist[1] < atlist[2]:
            out = [atlist[1]]
            if atlist[0] < atlist[2]:
                out.extend((atlist[0], atlist[2]))
            else:
                out.extend((atlist[2], atlist[0]))
        else:
            out = [atlist[2], atlist[1], atlist[0]]

    return out


def twosort_KvA(atlist):
    """
    Sorts a list of two items
    """
    if atlist[0] > atlist[1]:
        return [atlist[1], atlist[0]]
    else:
        return [atlist[0], atlist[1]]


# sorts a list of six items
def sixsort_KvA(atlist):
    """
    Sorts a list of six items
    """
    return union_KvA(threesort_KvA(atlist[:3]), threesort_KvA(atlist[3:]))
