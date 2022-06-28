# 3rd party lib imports
from numba import njit
import numpy as np

# own lib imports
from . import DataConversion as AIM_DC
from . import MathFunctions as AIM_MF
from . import PrintCommands as AIM_PC
from . import PhysicsFunctions as AIM_PF


def DoCoupling(FILES, RunPar, WS):
    """
    Called by CalcHam (AIM_PF), this function is the base of all coupling
    calculations for a single frame. First it prepares by doing some per-group
    calculations (like finding the relevant dipole moments). Then, it performs
    the coupling calculations themselves.
    """
    AIM_PC.vprint(
        4, ("Starting AIM_CF.PrepCoupling"),
        FILES.logfilename, RunPar)
    PrepCoupling(FILES, RunPar, WS)

    AIM_PC.vprint(
        4, ("Starting AIM_CF.CalcCoupling"),
        FILES.logfilename, RunPar)
    CalcCoupling(FILES, RunPar, WS)


def PrepCoupling(FILES, RunPar, WS):
    """
    Prepares the coupling. WS.coup_used is a dictionary with coupmap ID's as
    the keys, and a list of all groups that undergo that type of coupling at
    least once as the values. This way, only groups that undergo a certain
    coupling method are prepared for it (preparing for the wrong method could
    lead to all kinds of errors).
    """
    # first, prep the coupling
    for ID in WS.coup_used:
        AIM_PC.vprint(
            4, ("Preparing for the method with ID " + str(ID)),
            FILES.logfilename, RunPar)

        # The coupling maps that don't require prepping
        if ID in [0, 101, 102]:
            pass
        elif ID == 1:  # Generic TDC
            if RunPar.use_c_lib:
                setattr(WS, "TDCgen_r_c",
                        AIM_DC.ctype2d(WS.TDCgen_r, 'float32'))
                setattr(WS, "Dipoles_c",
                        AIM_DC.ctype2d(WS.Dipoles, 'float32'))
            else:
                pass
        elif ID == 2:  # TDCKrimm
            Prep_TDCKrimm(FILES, RunPar, WS)
        elif ID == 3:  # TDCTasumi
            Prep_TDCTasumi(FILES, RunPar, WS)
        elif ID == 4:  # TCC
            Prep_TCC(FILES, RunPar, WS)
        else:
            coupmap = RunPar.CouplingMaps[ID]
            coupmap.functions["prep_coupling"](
                FILES, RunPar, WS, coupmap)


def CalcCoupling(FILES, RunPar, WS):
    """
    Calculates all the couplings for a single frame. For each coupling,
    determines the method to use, and then applies it.
    """
    # This main loop is slow, as looping over the np array (by i and j) is
    # just incredibly slow. The only way to speed this up, is to not execute
    # this loop in python.
    # This is possible, though: with reason, couplings with nb/c support have
    # a low ID. Feed the coupling array to c (or nb, referred to as just c),
    # and let it check which function to apply. It calculates all the couplings
    # it can. The idea is then that python only has to fill the missing ones.
    # However, this still has a shortcoming: which entries are missing? Only
    # way to make that fast(er) is to create a sparse array of the i,j position
    # of all entries that have to be treated by python. This sparse array
    # can just be created once (it shouldn't change), or given by c to python.
    for OscGroupi in range(WS.res_desired_len):
        for OscGroupj in range(OscGroupi):
            CoupID = WS.coup_type_array[OscGroupi, OscGroupj]

            # this should just be 0. 1 = generic TDC (func still has to be
            # written and prepared!)
            if CoupID == 0:
                pass
            elif CoupID == 1:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcGenCoup(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            elif CoupID == 2:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcTDCKr(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            elif CoupID == 3:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcTDCTa(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            elif CoupID == 4:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcTCC(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            elif CoupID == 101:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcTasumi(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            elif CoupID == 102:
                WS.Hamiltonian[OscGroupj, OscGroupi] = CalcGLDP(
                    OscGroupi, OscGroupj, FILES, RunPar, WS)
            else:
                coupmap = RunPar.CouplingMaps[CoupID]
                WS.Hamiltonian[OscGroupj, OscGroupi] = (
                    coupmap.functions["calc_coupling"](
                        OscGroupi, OscGroupj, FILES, RunPar, WS, coupmap))


# ----


def Prep_TCC(FILES, RunPar, WS):
    """
    Prepares for calculating couplings using the TCC method.
    """
    alphaPro = WS.AllMaps.TCC.Pro_alpha
    alphaGen = WS.AllMaps.TCC.Gen_alpha
    if RunPar.use_c_lib:
        vPro = WS.AllMaps.TCC.Pro_v_c
        vGen = WS.AllMaps.TCC.Gen_v_c
        Prep_TCC_c(alphaPro, alphaGen, vPro, vGen, FILES.clib, WS)
    else:
        vPro = WS.AllMaps.TCC.Pro_v
        vGen = WS.AllMaps.TCC.Gen_v
        relevant_groups = WS.coup_used[4]
        WS.TCC_v = Prep_TCC_nb(
            alphaPro, alphaGen, vPro, vGen, WS.res_desired_len,
            relevant_groups,
            WS.AllOscGroups, WS.positions, WS.halfbox, WS.boxdims, WS.PrePro)


def Prep_TCC_c(alphaPro, alphaGen, vPro, vGen, clib, WS):
    """
    The c way of preparing for calculating couplings using TCC
    """
    clib.PrepTCC(
        WS.AllOscGroups_c, WS.PrePro_c, WS.TCC_v_c,
        WS.coup_used_c[4], len(WS.coup_used_c[4]), vPro,
        vGen, alphaPro, alphaGen, WS.positions_c, WS.halfbox_c, WS.boxdims_c
    )


@njit
def Prep_TCC_nb(alphaPro, alphaGen, vPro, vGen, res_desired_len,
                relevant_groups, AllAmGroups,
                WS_positions, WS_halfbox, WS_boxdims, WS_PrePro):
    """
    The numba way of preparing for calculating couplings using TCC
    """
    TCC_v = np.empty((res_desired_len, 18))
    for AmGroup in relevant_groups:
        C = AllAmGroups[AmGroup, 6]
        Ox = AllAmGroups[AmGroup, 7]
        N = AllAmGroups[AmGroup, 9]
        CO = AIM_MF.PBC_diff(WS_positions[Ox, :], WS_positions[C, :],
                             WS_halfbox, WS_boxdims)
        CO /= AIM_MF.vec3_len(CO)
        CN = AIM_MF.PBC_diff(WS_positions[N, :], WS_positions[C, :],
                             WS_halfbox, WS_boxdims)
        CN = AIM_MF.project(CO, CN)
        CN /= AIM_MF.vec3_len(CN)
        z = AIM_MF.crossprod(CO, CN)

        if WS_PrePro[AmGroup] == 0:
            for a in range(6):
                for c in range(3):
                    TCC_v[AmGroup, a*3 + c] = (
                        CO[c]*vGen[a, 0]
                        + CN[c]*vGen[a, 1]
                        + z[c]*vGen[a, 2])*alphaGen
        else:
            for a in range(6):
                for c in range(3):
                    TCC_v[AmGroup, a*3 + c] = (
                        CO[c]*vPro[a, 0]
                        + CN[c]*vPro[a, 1]
                        + z[c]*vPro[a, 2])*alphaPro
    return TCC_v


def Prep_TDCKrimm(FILES, RunPar, WS):
    """
    Prepares for calculating couplings using the Krimm transition dipole method
    """
    if RunPar.use_c_lib:
        Prep_TDCKrimm_c(FILES.clib, WS)
    else:
        relevant_groups = WS.coup_used[2]
        WS.TDCKr_r, WS.TDCKr_m = Prep_TDCKrimm_nb(
            WS.res_desired_len, relevant_groups, WS.AllOscGroups,
            WS.positions, WS.halfbox, WS.boxdims)


def Prep_TDCKrimm_c(clib, WS):
    """
    The c way of preparing for calculating couplings for the Krimm transition
    dipole method.
    """
    clib.PrepTDCKrimm(
        WS.AllOscGroups_c, WS.TDCKr_r_c, WS.TDCKr_m_c,
        WS.coup_used_c[2], len(WS.coup_used_c[2]), WS.positions_c,
        WS.halfbox_c, WS.boxdims_c)


@njit
def Prep_TDCKrimm_nb(res_desired_len, relevant_groups, AllAmGroups,
                     WS_positions, halfbox, boxdims):
    """
    The numba way of preparing for calculating couplings for the Krimm
    transition dipole method.
    """
    displace = 0.868
    tanangle = -0.36397023426

    TDCKr_r = np.empty((res_desired_len, 3))
    TDCKr_m = np.empty((res_desired_len, 3))

    for AmGroup in relevant_groups:
        C = AllAmGroups[AmGroup, 6]
        Ox = AllAmGroups[AmGroup, 7]
        N = AllAmGroups[AmGroup, 9]

        CO = AIM_MF.PBC_diff(WS_positions[Ox, :], WS_positions[C, :],
                             halfbox, boxdims)

        dispibond = displace / AIM_MF.vec3_len(CO)
        TDCKr_r[AmGroup, :] = WS_positions[C, :] + dispibond * CO
        CO /= AIM_MF.vec3_len(CO)

        CN = AIM_MF.PBC_diff(WS_positions[N, :], WS_positions[C, :],
                             halfbox, boxdims)
        CN = AIM_MF.project(CO, CN)
        CN /= AIM_MF.vec3_len(CN)

        bond = tanangle
        m = CO + CN*bond
        m /= AIM_MF.vec3_len(m)
        TDCKr_m[AmGroup, :] = m

    return TDCKr_r, TDCKr_m


def Prep_TDCTasumi(FILES, RunPar, WS):
    """
    Prepares for calculating couplings using the Tasumi transition dipole
    method.
    """
    if RunPar.use_c_lib:
        Prep_TDCTasumi_c(FILES.clib, WS)
    else:
        relevant_groups = WS.coup_used[3]
        WS.TDCTa_r, WS.TDCTa_m = Prep_TDCTasumi_nb(
            WS.res_desired_len, relevant_groups, WS.AllOscGroups, WS.positions,
            WS.halfbox, WS.boxdims)


def Prep_TDCTasumi_c(clib, WS):
    """
    The c way of preparing for calculating couplings for the Tasumi transition
    dipole method.
    """
    clib.PrepTDCTasumi(
        WS.AllOscGroups_c, WS.TDCTa_r_c, WS.TDCTa_m_c,
        WS.coup_used_c[3], len(WS.coup_used_c[3]), WS.positions_c,
        WS.halfbox_c, WS.boxdims_c)


@njit
def Prep_TDCTasumi_nb(res_desired_len, relevant_groups, AllAmGroups,
                      WS_positions, halfbox, boxdims):
    """
    The numba way of preparing for calculating couplings for the Tasumi
    transition dipole method.
    """
    itheta = 5.6712818196
    # = 1/tan(0.17632698 = 10 deg expressed in rad)

    # magnitude = 2.73 * np.sqrt(51.43 / 5034)
    magnitude = 0.276

    TDCTa_r = np.empty((res_desired_len, 3))
    TDCTa_m = np.empty((res_desired_len, 3))
    for AmGroup in relevant_groups:
        C = AllAmGroups[AmGroup, 6]
        Ox = AllAmGroups[AmGroup, 7]
        N = AllAmGroups[AmGroup, 9]
        CO = AIM_MF.PBC_diff(WS_positions[Ox, :], WS_positions[C, :],
                             halfbox, boxdims)
        CN = AIM_MF.PBC_diff(WS_positions[N, :], WS_positions[C, :],
                             halfbox, boxdims)
        CO /= AIM_MF.vec3_len(CO)
        CN /= AIM_MF.vec3_len(CN)
        dri = 0.665*CO + 0.258*CN
        TDCTa_r[AmGroup, :] = WS_positions[C, :]+dri

        COvecDri = AIM_MF.dotprod(CO, dri)
        prefac = COvecDri + np.sqrt(AIM_MF.dotprod(dri, dri)
                                    - COvecDri*COvecDri)*itheta
        m = dri - prefac*CO
        m /= AIM_MF.vec3_len(m)
        TDCTa_m[AmGroup, :] = m*magnitude
    return TDCTa_r, TDCTa_m


# ----


def CalcGenCoup(OscGroupi, OscGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the general
    dipole-dipole coupling. Positions should have units of angstrom, dipole
    moments should have units of debye
    """
    if RunPar.use_c_lib:
        J = FILES.clib.CalcGenCoup(
            OscGroupi, OscGroupj, WS.TDCgen_r_c, WS.Dipoles_c,
            WS.halfbox_c, WS.boxdims_c)
    else:
        J = CalcGenCoup_nb(
            OscGroupi, OscGroupj, WS.TDCgen_r, WS.Dipoles, WS.halfbox,
            WS.boxdims)

    J *= RunPar.Scale_LRCoupling
    return J


@njit
def CalcGenCoup_nb(OscGroupi, OscGroupj, TDCgen_r, TDCgen_m, halfbox, boxdims):
    """
    The numba way of calculating the coupling between two given groups using
    the generic dipole-dipole coupling approximation
    """
    # Used constants:
    # Cm = (1/3.33564) * 10^30 D  (Coulomb meter in Debye)
    # m = 10^10 ang (meter in angstrom)
    # J = 1/hc = (1/1.98644586) * 10^25 1/m
    # => J = 5.03411656 * 10^22 1/cm (joule in wavenumbers)
    # eps_0 = 8.8541878128 F/m = 8.8541878128 C^2/Jm (coulomb squared per
    # joule meter)

    # derived value:
    # 4piEinv = 1/(4 * pi * eps_0) Jm/C^2
    # Gives 5034.11656 cm^-1 * ang*3 Deb^-2

    fourPiEps_inv = 5034
    d = AIM_MF.PBC_diff(TDCgen_r[OscGroupi, :], TDCgen_r[OscGroupj, :],
                        halfbox, boxdims)
    ir2 = 1/AIM_MF.dotprod(d, d)
    ir = np.sqrt(ir2)
    ir3 = ir*ir2
    ir5 = ir3*ir2

    J = fourPiEps_inv * (
        AIM_MF.dotprod(TDCgen_m[OscGroupi], TDCgen_m[OscGroupj]) * ir3
        - 3.0 * AIM_MF.dotprod(TDCgen_m[OscGroupi], d)
        * AIM_MF.dotprod(TDCgen_m[OscGroupj], d) * ir5)

    return J


def CalcTCC(AmGroupi, AmGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the TCC method.
    """
    if RunPar.use_c_lib:
        J = CalcTCC_c(AmGroupi, AmGroupj, FILES.clib, WS)
    else:
        J = CalcTCC_nb(AmGroupi, AmGroupj, WS.positions, WS.AllOscGroups,
                       WS.halfbox, WS.boxdims, WS.PrePro, WS.AllMaps.TCC.Gen_q,
                       WS.AllMaps.TCC.Pro_q, WS.AllMaps.TCC.Gen_dq,
                       WS.AllMaps.TCC.Pro_dq, WS.TCC_v,
                       WS.AllMaps.TCC.fourPiEps)

    J *= RunPar.Scale_LRCoupling

    return J


def CalcTCC_c(AmGroupi, AmGroupj, clib, WS):
    """
    The c way of calculating the coupling between two given groups using the
    TCC method.
    """
    J = clib.CalcTCC(AmGroupi, AmGroupj, WS.AllOscGroups_c, WS.PrePro_c,
                     WS.TCC_v_c, WS.AllMaps.TCC.Gen_q_c,
                     WS.AllMaps.TCC.Pro_q_c, WS.AllMaps.TCC.Gen_dq_c,
                     WS.AllMaps.TCC.Pro_dq_c, WS.AllMaps.TCC.fourPiEps,
                     WS.positions_c, WS.halfbox_c, WS.boxdims_c)
    return J


@njit
def CalcTCC_nb(AmGroupi, AmGroupj, WS_positions, AllAmGroups, halfbox, boxdims,
               PrePro, qGen, qPro, dqGen, dqPro, TCC_v, TCC_4PiEps):
    """
    The numba way of calculating the coupling between two given groups using
    the TCC method.
    """
    J = 0

    tyi = PrePro[AmGroupi]
    tyj = PrePro[AmGroupj]

    qi = tyi * qPro + (1-tyi) * qGen
    dqi = tyi * dqPro + (1-tyi) * dqGen
    qj = tyj * qPro + (1-tyj) * qGen
    dqj = tyj * dqPro + (1-tyi) * dqGen

    for a in range(6):
        va = TCC_v[AmGroupi, a*3:(a+1)*3]
        for b in range(6):
            d = AIM_MF.PBC_diff(
                WS_positions[AllAmGroups[AmGroupi, 6+a], :],
                WS_positions[AllAmGroups[AmGroupj, 6+b], :],
                halfbox, boxdims)
            r2 = AIM_MF.dotprod(d, d)
            if r2 < 0.01:
                r2 = 1

            ir2 = 1/r2
            ir = np.sqrt(ir2)
            ir3 = ir * ir2
            ir5 = ir3 * ir2

            vb = TCC_v[AmGroupj, b*3:(b+1)*3]
            J -= (3 * ir5
                  * qi[a]
                  * qj[b]
                  * AIM_MF.dotprod(vb, d)
                  * AIM_MF.dotprod(va, d))
            J -= ir3 * (dqi[a] * qj[b] * AIM_MF.dotprod(vb, d)
                        + qi[a] * dqj[b] * AIM_MF.dotprod(va, d)
                        - AIM_MF.dotprod(va, vb) * qi[a] * qj[b])
            J += ir * dqi[a] * dqj[b]
    J *= TCC_4PiEps

    return J


def CalcTDCKr(AmGroupi, AmGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the Krimm
    transition dipole method.
    """
    if RunPar.use_c_lib:
        J = FILES.clib.CalcTDCKrimm(
            AmGroupi, AmGroupj, WS.TDCKr_r_c, WS.TDCKr_m_c,
            WS.halfbox_c, WS.boxdims_c)
    else:
        J = CalcTDCKr_nb(
            AmGroupi, AmGroupj, WS.TDCKr_r, WS.TDCKr_m,
            WS.halfbox, WS.boxdims)

    J *= RunPar.Scale_LRCoupling

    return J


@njit
def CalcTDCKr_nb(AmGroupi, AmGroupj, TDCKr_r, TDCKr_m, halfbox, boxdims):
    """
    The numba way of calculating the coupling between two given groups using
    the Krimm transition dipole method.
    """
    fourPiEps = 580
    d = AIM_MF.PBC_diff(TDCKr_r[AmGroupi, :], TDCKr_r[AmGroupj, :],
                        halfbox, boxdims)
    ir2 = 1/AIM_MF.dotprod(d, d)
    ir = np.sqrt(ir2)
    ir3 = ir*ir2
    ir5 = ir3*ir2
    J = fourPiEps * (
        AIM_MF.dotprod(TDCKr_m[AmGroupi], TDCKr_m[AmGroupj]) * ir3
        - 3.0 * AIM_MF.dotprod(TDCKr_m[AmGroupi], d)
        * AIM_MF.dotprod(TDCKr_m[AmGroupj], d) * ir5)

    return J


def CalcTDCTa(AmGroupi, AmGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the Tasumi
    transition dipole method.
    """
    if RunPar.use_c_lib:
        J = FILES.clib.CalcTDCTasumi(
            AmGroupi, AmGroupj, WS.TDCTa_r_c, WS.TDCTa_m_c,
            WS.halfbox_c, WS.boxdims_c)
    else:
        J = CalcTDCTa_nb(
            AmGroupi, AmGroupj, WS.TDCTa_r, WS.TDCTa_m,
            WS.halfbox, WS.boxdims)
    J *= RunPar.Scale_LRCoupling
    return J


@njit
def CalcTDCTa_nb(AmGroupi, AmGroupj, TDCTa_r, TDCTa_m, halfbox, boxdims):
    """
    The numba way of calculating the coupling between two given groups using
    the Tasumi transition dipole method.
    """
    # AoverEps = 51.43
    AoverEps = 5034  # see CalcGenCoup

    d = AIM_MF.PBC_diff(TDCTa_r[AmGroupi, :], TDCTa_r[AmGroupj, :],
                        halfbox, boxdims)
    ir2 = 1/AIM_MF.dotprod(d, d)
    ir = np.sqrt(ir2)
    ir3 = ir * ir2
    ir5 = ir3 * ir2
    J = AoverEps * (
        AIM_MF.dotprod(TDCTa_m[AmGroupi], TDCTa_m[AmGroupj]) * ir3
        - 3.0 * AIM_MF.dotprod(TDCTa_m[AmGroupi], d)
        * AIM_MF.dotprod(TDCTa_m[AmGroupj], d) * ir5)

    return J


def CalcTasumi(AmGroupi, AmGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the Tasumi method.
    This method is specialized for through-bond coupling of amides in the
    protein backbone.
    """
    if WS.AllOscGroups[AmGroupi, 2] == WS.AllOscGroups[AmGroupj, 8]:
        bond = WS.AllOscGroups[AmGroupi, :12]
    elif WS.AllOscGroups[AmGroupi, 12] == WS.AllOscGroups[AmGroupj, 8]:
        bond = WS.AllOscGroups[AmGroupi, 6:]

    phi_ang = AIM_MF.dihedral(
        WS.positions[bond[0], :], WS.positions[bond[3], :],
        WS.positions[bond[5], :], WS.positions[bond[6], :],
        WS.halfbox, WS.boxdims
        ) * (180/3.1416)
    psi_ang = AIM_MF.dihedral(
        WS.positions[bond[3], :], WS.positions[bond[5], :],
        WS.positions[bond[6], :], WS.positions[bond[9], :],
        WS.halfbox, WS.boxdims
        ) * (180/3.1416)

    phi_N = int((phi_ang+180)//30)
    psi_N = int((psi_ang+180)//30)

    x1l = phi_N * 30 - 180
    x2l = psi_N * 30 - 180

    y1 = WS.AllMaps.Tasumi.map[phi_N, psi_N]
    y2 = WS.AllMaps.Tasumi.map[phi_N+1, psi_N]
    y3 = WS.AllMaps.Tasumi.map[phi_N+1, psi_N+1]
    y4 = WS.AllMaps.Tasumi.map[phi_N, psi_N+1]

    u = (psi_ang - x2l)/30
    t = (phi_ang - x1l)/30

    J = (1-u)*(1-t)*y1 + (1-u)*t*y2 + u*t*y3 + u*(1-t)*y4
    return J


def CalcGLDP(AmGroupi, AmGroupj, FILES, RunPar, WS):
    """
    Calculates the coupling for a given pair of groups using the GLDP method.
    This method is specialized for through-bond coupling of amides in the
    protein backbone.
    """
    if WS.AllOscGroups[AmGroupi, 2] == WS.AllOscGroups[AmGroupj, 8]:
        bond = 0
    elif WS.AllOscGroups[AmGroupi, 12] == WS.AllOscGroups[AmGroupj, 8]:
        bond = 1

    mapname = AIM_PF.cis_or_trans(
        WS.AllOscGroups[AmGroupi, :], bond+1, RunPar, WS)
    NNmap = getattr(WS.AllMaps.NN, "Coupling" + mapname)

    J = AIM_PF.calcdelta(
        WS.AllOscGroups[AmGroupi, :], NNmap, bond*6, 0,
        FILES.logfilename, RunPar, WS)[0]

    return J
