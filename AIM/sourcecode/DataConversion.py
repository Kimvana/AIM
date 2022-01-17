# 3rd party lib imports
from numba import njit
import numpy as np


def ctype2d(array, type):
    """
    Takes a 2D numpy array and converts it to a c-friendly version of the
    desired type. These arrays require a change of dimensions (done by
    np.ravel).
    """
    array = array.astype(type)
    arr = np.ravel(array)
    out = np.ctypeslib.as_ctypes(arr)
    return out


def ctypelist(array, type):
    """
    Takes a python list and converts it to a c-friendly array of the desired
    type.
    """
    arr = np.array(array).astype(type)
    out = np.ctypeslib.as_ctypes(arr)
    return out


def timestring(time):
    """
    takes a time (float) in seconds, and converts this to a day-hour:min:sec
    notation, as a string.
    """
    time_sec = str(round(time % 60, 3))
    time_sec_check = str(int(float(time_sec)))
    time = int(time//60)
    time_min = str(time % 60)
    time = time//60
    time_hour = str(time % 24)
    time_day = str(time//24)

    time_sec = (2-len(time_sec_check))*"0" + time_sec
    time_min = (2-len(time_min))*"0" + time_min
    time_hour = (2-len(time_hour))*"0" + time_hour
    time_day = (3-len(time_day))*" " + time_day

    string = time_day + "-" + time_hour + ":" + time_min + ":" + time_sec
    return string


def ix_builder_c(inrange_reslist_c, inrange_reslist_len, clib, WS):
    """
    Wrapper for the c code for converting a list of residue numbers to a list
    of the numbers of the atoms in those residues.
    """
    outlen = clib.ix_builder(
        inrange_reslist_c, inrange_reslist_len, WS.residues.start_c,
        WS.residues.fin_c, WS.ix.inrange_long_c)
    inrange_ix_c = np.ctypeslib.as_array(WS.ix.inrange_long_c)
    return inrange_ix_c[:outlen]


@njit
def posfinder(inrange_ix, WS_positions, WS_charges):
    """
    Slices the general WS position and charge arrays to only include the
    required atoms (the atoms stored in inrange_ix)
    """
    inrange_natoms = len(inrange_ix)
    inrange_pos = np.zeros((inrange_natoms, 3))
    inrange_char = np.zeros((inrange_natoms))
    for index, atnum in enumerate(inrange_ix):
        inrange_pos[index, :] = WS_positions[atnum, :]
        inrange_char[index] = WS_charges[atnum]
    return inrange_pos, inrange_char
