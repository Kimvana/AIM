# 3rd party lib imports
from numba import njit
import numpy as np


@njit
def PBC_diff(vect1, vect2, halfbox, boxdims):
    """
    Calculates the distance vector between vect1 and vect2, taking into account
    the periodic nature of the system. diff = vect1 - vect2 (diff = the vector
    pointing from vect2 to vect1).
    """
    diff = np.copy(vect1)
    for i in range(3):
        diff[i] -= vect2[i]
        if diff[i] > halfbox[i]:
            diff[i] -= boxdims[i]
        elif diff[i] < -1*halfbox[i]:
            diff[i] += boxdims[i]
    return diff


@njit
def crossprod(vect1, vect2):
    """
    Calculates the cross product between two vectors of size 3. This is
    faster than the dedicated np method, as there are no checks for the
    correctness of the provided vectors.
    """
    vect3 = np.copy(vect1)
    vect3[0] = vect1[1]*vect2[2]-vect1[2]*vect2[1]
    vect3[1] = vect1[2]*vect2[0]-vect1[0]*vect2[2]
    vect3[2] = vect1[0]*vect2[1]-vect1[1]*vect2[0]

    return vect3


@njit
def dotprod(vect1, vect2):
    """
    Calculates the dot product between two vectors of size 3. This is faster
    than the dedicated np method, as there are no checks for the correctness
    of the provided vectors.
    """
    return vect1[0]*vect2[0] + vect1[1]*vect2[1] + vect1[2]*vect2[2]


@njit
def vec3_len(vect):  # replacement for np.linalg.norm
    """
    Calculates the norm (length) of a vector of size 3. This is faster than the
    dedicated np method, as there are no checks for the correctness of the
    provided vectors.
    """
    return np.sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2])


@njit
def project(a, b):
    """
    Calculates the part of vector b that is orthogonal to the vector a (i.e.
    it subtracts from b the part that is along a, and returns the result).
    """
    inprod = dotprod(a, b)/dotprod(a, a)
    b -= inprod*a
    return b


def tensorprod(vect1):
    """
    Calculates the tensor product of two identical vectors of size 3.
    Returns the matrix as a vector of length 6 (xx xy xz yy yz zz)
    """
    tensor = np.zeros((6), dtype='float32')
    tensor[:3] = vect1[0]*vect1
    tensor[3:5] = vect1[1]*vect1[1:]
    tensor[5] = vect1[2]*vect1[2]

    return tensor


@njit
def dihedral(p0, p1, p2, p3, halfbox, boxdims):
    """
    Calculates the dihedral angle between the supplied points. The order of the
    points matters: the calculated angle is the following:
    When looking at these 4 points such that p2 lies behind p1 (or other way
    around?), we're concerned with the apparent angle p0, p1, p3 (which is
    identical to the apparent angle p0, p2, p3).
    """
    b0 = PBC_diff(p0, p1, halfbox, boxdims)
    b1 = PBC_diff(p2, p1, halfbox, boxdims)
    b2 = PBC_diff(p3, p2, halfbox, boxdims)

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= vec3_len(b1)

    # = projection of b0 onto plane perpendicular to b1 (= b0 minus component
    # that aligns with b1)
    # In principle, v = project(b1, b0) would give the same result. However,
    # project is more expensive, as it divides by dot(b1, b1). This is
    # basically normalizing, which is here done prior already.
    v = b0 - dotprod(b0, b1)*b1

    # = projection of b2 onto plane perpendicular to b1 (= b2 minus component
    # that aligns with b1)
    w = b2 - dotprod(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x

    # np.arctan2: computes angle between the vector pointing to (x, y) and the
    # vector (1, 0) (x axis). arg1 = y, arg2 = x.
    x = dotprod(v, w)  # how much v and w align
    y = dotprod(crossprod(b1, v), w)  # cross rotates v 90 degrees.
    return np.arctan2(y, x)


@njit
def EG_rotation(E, G, U, relevant_j):
    """
    rotates the fields and gradients to the correct direction specified by the
    rotation matrix U.

    E and G are expected to have the default format: like they are supplied
    when generated. This means both have 6 rows (one for each atom). E has
    three columns, corresponding to the x, y, z directions resp. G has six
    columns, corresponding to the xx, yy, zz, xy, xz, yz directions resp.

    U is the (3 by 3) matrix to convert from the box-relative coordinates to
    the oscillator-relative coordinates such that E_osc = U (dot) E_box.
    It's three rows are the new, oscillator-relative x, y and z vectors,
    each expressed in box-relative coordinates.
    Imagine a rotation where the new X axis = old Z, new Y = old X, new Z = old
    Y.      (0 0 1)
        U = (1 0 0)
            (0 1 0)
    """

    # old version:
    # for j in relevant_j:
    #     E[j, :] = np.dot(rot_mat, E[j, :])

    #     Gsq = np.diag(G[j, :3])
    #     Gsq[np.triu_indices(3, k=1)] = G[j, 3:]
    #     Gsq[np.tril_indices(3, k=-1)] = G[j, 3:]
    #     temp = np.matmul(
    #         np.matmul(rot_mat, Gsq),
    #         rot_mat.transpose()
    #     )
    #     G[j, :3] = np.diag(temp)
    #     G[j, 3:] = temp[np.triu_indices(3, k=1)]

    for j in relevant_j:
        E[j, :] = np.dot(U, E[j, :])

        Gsq = np.diag(G[j, :3])
        Gsq[0, 1:] = G[j, 3:5]
        Gsq[1, 2] = G[j, 5]
        Gsq[1:, 0] = G[j, 3:5]
        Gsq[2, 1] = G[j, 5]
        temp = U @ Gsq @ U.transpose()
        G[j, :3] = np.diag(temp)
        G[j, 3:5] = temp[0, 1:]
        G[j, 5] = temp[1, 2]

    return E, G


@njit
def apply_map(P, E, G, omega, Pmap, Emap, Gmap, relevant_j):
    """
    A function for applying 'standard' maps. Allows for any number of atoms,
    but assumes a 'standard' P, E, and G map. This kind of structure should be
    possible for any non-symmetric system.
    """
    returnval = omega

    returnval += np.sum(np.multiply(P[relevant_j], Pmap[relevant_j]))
    returnval += np.sum(np.multiply(E[relevant_j, :], Emap[relevant_j, :]))
    returnval += np.sum(np.multiply(G[relevant_j, :], Gmap[relevant_j, :]))

    return returnval


@njit
def apply_map_1D(P, Pmap, relevant_j):
    """
    A more general function for applying maps. This is for 1D maps, like for
    potential.
    """
    return np.sum(np.multiply(P[relevant_j], Pmap[relevant_j]))


@njit
def apply_map_2D(E, Emap, relevant_j):
    """
    A more general function for applying maps. This is for 2D maps, like for
    field/gradient.
    """
    return np.sum(np.multiply(E[relevant_j, :], Emap[relevant_j, :]))
