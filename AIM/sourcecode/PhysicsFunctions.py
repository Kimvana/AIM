# 3rd party lib imports:
from numba import njit
import numpy as np


# my lib imports:
from . import CouplingFunctions as AIM_CF
from . import DataConversion as AIM_DC
from . import MathFunctions as AIM_MF
from . import SetOperations as AIM_SO
from . import PrintCommands as AIM_PC


def calc_COM(choice, FILES, RunPar, WS):
    """
    Calculates the centre of mass (COM) of residues. The input argument
    'choice' determines which: if choice == "Sel", the COM of each residue in
    the MD system is calculated (this is the slow one). If the user opts to
    use NSA, this is one of the two steps where it matters: Every n frames,
    NSA still calls this function using "Sel" (then, it updates the master
    list). All the other, it calls it using "NSA". In those cases, only the COM
    for residues that are within a range of NSA_spheresize is calculated. This
    is what determines the choice arrays (choice_c, choice_len, choice_ar).
    These arrays contain the residues we're interested in.

    Then, there are two ways to calculate the COM. One does take into account
    the perodic nature of the simulation box, the other doesn't. Of course, the
    first is physically accurate. But the second one replicates the results of
    the original AIM script. Therefore, the faulty one (dubbed 'dumb') is
    used if the user requested replication of the original behaviour.
    """
    if choice == "NSA":
        if RunPar.use_c_lib:
            choice_c = WS.all_inrange_res_c
            choice_len = WS.all_inrange_nres
        else:
            choice_ar = WS.all_inrange_res
    elif choice == "Sel":
        if RunPar.use_c_lib:
            choice_c = WS.residues.resnums_COM_c
            choice_len = WS.residues.resnums_COM_len
        else:
            choice_ar = WS.residues.resnums_COM

    if RunPar.replicate_orig_AIM:
        if RunPar.use_c_lib:
            Res_COM = calc_COM_dumb_c(WS, FILES.clib, choice_c, choice_len)
        else:
            Res_COM = calc_COM_dumb_nb(
                choice_ar, WS.nres, WS.residues.start, WS.residues.fin,
                WS.positions, WS.masses)
    else:
        if RunPar.use_c_lib:
            Res_COM = calc_COM_c(WS, FILES.clib, choice_c, choice_len)
        else:
            Res_COM = calc_COM_nb(
                choice_ar, WS.nres, WS.residues.start, WS.residues.fin,
                WS.positions, WS.masses, WS.halfbox, WS.boxdims)

    return Res_COM


@njit
def calc_COM_nb(COM_resnums, nres, resstart, resfin, positions, masses,
                halfbox, boxdims):
    """
    The numba way to correctly calculate the COM's.
    """
    res_array = np.zeros((nres, 3))
    for resnum in COM_resnums:
        start = resstart[resnum]
        fin = resfin[resnum]

        # sum up mass and positions of all atoms (taking into account the PBC)
        refpos = positions[start, :]
        cumpos = np.zeros((3))
        cummass = masses[start]
        for atom in range(start+1, fin+1):
            diff = AIM_MF.PBC_diff(positions[atom, :], refpos,
                                   halfbox, boxdims)
            cumpos += (diff*masses[atom])
            cummass += masses[atom]

        # divide the two
        res_array[resnum, :] = (cumpos/cummass) + refpos

        # if too far out -> place it back!
        for i in range(3):
            if res_array[resnum, i] > halfbox[i]:
                res_array[resnum, i] -= boxdims[i]
            elif res_array[resnum, i] < (-1*halfbox[i]):
                res_array[resnum, i] += boxdims[i]
    return res_array


def calc_COM_c(WS, clib, resnums_c, resnums_len):
    """
    The c-way to correctly calculate the COM's.
    """
    # call the c-script, calculate resnums
    clib.res_finder(resnums_c, resnums_len, WS.residues.start_c,
                    WS.residues.fin_c, WS.positions_c, WS.masses_c,
                    WS.res_COM_c, WS.halfbox_c, WS.boxdims_c)

    # convert COM_array back to more easily readable format!
    res_COM = np.ctypeslib.as_array(WS.res_COM_c)
    COM_array = np.empty((WS.nres, 3))
    for i in range(WS.nres):
        COM_array[i, :] = res_COM[i*3:(i+1)*3]
    return COM_array


@njit
def calc_COM_dumb_nb(COM_resnums, nres, resstart, resfin, positions, masses):
    """
    the numba way to calculate the COM's without taking into account the pbc.
    """
    res_array = np.zeros((nres, 3))
    for resnum in COM_resnums:
        start = resstart[resnum]
        fin = resfin[resnum]

        # sum up mass and positions of all atoms (taking into account the PBC)
        cumpos = np.zeros((3))
        cummass = 0
        for atom in range(start, fin+1):
            diff = positions[atom, :]
            cumpos += (diff*masses[atom])
            cummass += masses[atom]

        # divide the two
        res_array[resnum, :] = (cumpos/cummass)

    return res_array


def calc_COM_dumb_c(WS, clib, resnums_c, resnums_len):
    """
    the c way to calculate the COM without taking into account the pbc.
    """
    # call the c-script, calculate resnums
    clib.res_finder_dumb(resnums_c, resnums_len, WS.residues.start_c,
                         WS.residues.fin_c, WS.positions_c, WS.masses_c,
                         WS.res_COM_c, WS.halfbox_c, WS.boxdims_c)

    # convert COM_array back to more easily readable format!
    res_COM = np.ctypeslib.as_array(WS.res_COM_c)
    COM_array = np.empty((WS.nres, 3))
    for i in range(WS.nres):
        COM_array[i, :] = res_COM[i*3:(i+1)*3]
    return COM_array


# ----------


def AGsorterNSA(refind, FILES, RunPar, WS):
    """
    Called by Universe.COM_update (AIM_UM), but only when using NSA. Takes all
    residues in the system, and sees which of those are within range of the
    residue with resum == refind (=reference index).
    """
    if RunPar.use_c_lib:
        returnobj = AGsorter_c(refind, RunPar.NSA_spheresize_c,
                               WS.residues.resnums_influencers_c,
                               WS.residues.resnums_influencers_len,
                               FILES.clib, WS)
    else:
        returnobj = AGsorter_nb(refind, WS.residues.resnums_influencers,
                                WS.Res_COM, WS.halfbox, WS.boxdims,
                                WS.residues.start, WS.residues.fin,
                                RunPar.NSA_spheresize)[1]
        returnobj = np.array(returnobj)
    return returnobj


def AGsorterEstat(refind, FILES, RunPar, WS):
    """
    Called by CalcDiag (AIM_PF). Takes a longer list of residues, and checks
    which of these are within range = SphereSize of the residue with index
    refind. If NSA is used, this longer list is the list found by AGsorterNSA,
    otherwise it's all the residues in the MD system.
    """
    if RunPar.use_c_lib:
        if RunPar.NSA_toggle:
            inrange_reslist = AGsorter_c(refind, RunPar.SphereSize_c,
                                         WS.NSA_resdict[refind],
                                         WS.NSA_lendict[refind],
                                         FILES.clib, WS)
        else:
            inrange_reslist = AGsorter_c(refind, RunPar.SphereSize_c,
                                         WS.residues.resnums_influencers_c,
                                         WS.residues.resnums_influencers_len,
                                         FILES.clib, WS)

        # AGsorter_c itself returns a list of residue numbers, not atom
        # numbers (is faster). ix_builder_c converts the residue number
        # list into an atom number one.
        inrange_ix = AIM_DC.ix_builder_c(
            np.ctypeslib.as_ctypes(inrange_reslist),
            np.int32(len(inrange_reslist)), FILES.clib, WS)
        inrange_ix = np.ctypeslib.as_ctypes(inrange_ix)
    else:
        if RunPar.NSA_toggle:
            inrange_ix, inrange_reslist = AGsorter_nb(
                refind, WS.NSA_resdict[refind], WS.Res_COM, WS.halfbox,
                WS.boxdims, WS.residues.start, WS.residues.fin,
                RunPar.SphereSize)
        else:
            inrange_ix, inrange_reslist = AGsorter_nb(
                refind, WS.residues.resnums_influencers, WS.Res_COM,
                WS.halfbox, WS.boxdims, WS.residues.start, WS.residues.fin,
                RunPar.SphereSize)

    return inrange_ix


def AGsorter_c(refind, desrange, influencers_residues_resnums_c,
               sel_res_resnums_len, clib, WS):
    """
    The c-way of calculating which residues are within range of the reference
    residue.
    """
    out = clib.AGsorter(refind, influencers_residues_resnums_c,
                        sel_res_resnums_len, WS.res_COM_c,
                        WS.residues.inrange_c, desrange, WS.nres,
                        WS.halfbox_c, WS.boxdims_c)

    # convert to be able to cut the array.
    inrange_res = np.ctypeslib.as_array(WS.residues.inrange_c)
    inrange_res_copy = np.copy(inrange_res)
    return inrange_res_copy[:out]


@njit
def AGsorter_nb(refind, influencers_residues_resnums, WS_res_COM, WS_halfbox,
                WS_boxdims, WS_residues_start, WS_residues_fin, desrange):
    """
    The numba way of calculating which residues are within range of the
    reference residue.
    """
    Am_COM = WS_res_COM[refind, :]
    inside_atnums, inside_residues_resnums = [], []
    desrange2 = desrange*desrange

    # for each of the residues, look at distance to reference position. If too
    # far, do not include!
    for resnum in influencers_residues_resnums:
        COM = WS_res_COM[resnum, :]
        diff = AIM_MF.PBC_diff(Am_COM, COM, WS_halfbox, WS_boxdims)
        dist2 = AIM_MF.dotprod(diff, diff)
        if dist2 < desrange2:
            inside_residues_resnums.append(resnum)
            inside_atnums.extend(
                range(WS_residues_start[resnum], WS_residues_fin[resnum]+1))
    return inside_atnums, inside_residues_resnums


# ---------


def CalcHam(FILES, RunPar, WS):
    """
    Called by Universe.RunCalc (AIM_UM). Calculates all the entries to the
    Hamiltonian. The diagonal elements require information on the electric
    field, just like the dipoles do. Therefore, the dipoles are also calculated
    here.
    """
    # writes diagonal entries directly into hamiltonian stored in WS, and
    # Dipoles into Dipoles stored in WS
    if any(
        outtype in RunPar.output_type for outtype in [
            "Ham", "Dip", "Ram"
        ]
    ):
        AIM_PC.vprint(
            4, ("Starting AIM_PF.CalcDiag"),
            FILES.logfilename, RunPar)
        CalcDiag(FILES, RunPar, WS)

    if "Ham" in RunPar.output_type:
        AIM_PC.vprint(
            4, ("Starting AIM_PF.DoCoupling"),
            FILES.logfilename, RunPar)
        AIM_CF.DoCoupling(FILES, RunPar, WS)

    if "Pos" in RunPar.output_type:
        AIM_PC.vprint(
            4, ("Starting AIM_PF.FindAtomPos"),
            FILES.logfilename, RunPar)

        FindAtomPos(WS)


def CalcDiag(FILES, RunPar, WS):
    """
    This function calculates the diagonal elements (oscillation frequencies) of
    the hamiltonian. As this requires the electric field properties, just like
    the dipoles, the dipoles are also calculated here
    """
    function_warning = 0
    for AmGroup in range(WS.res_desired_len):
        Am_C = WS.AllOscGroups[AmGroup, 6]
        resnum = WS.resnums[Am_C]

        # We want to calculate the Estatic effects only for the atoms that
        # fall within the specified SphereSize range. AGsorter returns only
        # the atom numbers that fall within the range out of the provided list.
        # In the end, obtain inrange_ix (datatype depends on use_c_lib), which
        # contains the atom numbers of all the atoms that are close enough for
        # Estatic interactions. If NSA_toggle, there is a shortlist of atom
        # numbers that is faster to loop over than the full list used
        # otherwise!
        inrange_ix = AGsorterEstat(resnum, FILES, RunPar, WS)

        P, E, G = CalcEstat(inrange_ix, AmGroup, FILES, RunPar, WS)

        # regardless of map_choice, if it is a sidechain group, we need the
        # sidechain map. map_reader made sure that we have the correct map
        # present. Now, based on the Estatic effects, calculate the diagonal
        # elements of the hamiltonian
        rot_mat = CalcHamDiagRaw(AmGroup, P, E, G, FILES, RunPar, WS)

        # nnmap directly alters the hamiltonian.
        if (
            "Ham" in RunPar.output_type
            and WS.OscID[AmGroup] == 0
            and RunPar.TreatNN
            and RunPar.map_choice != "Tokmakoff"
        ):
            function_warning = nnmap(AmGroup, function_warning,
                                     FILES.logfilename, RunPar, WS)

        # this is only relevant to backbone amide groups. The rest is solved
        # by the .map files!
        # why? Jansen dipoles should always return the same dipoles, but it's
        # dependent on the field created... This makes them independent of
        # map_choice.
        if WS.OscID[AmGroup] == 0:
            P, E, G = UpdateEstat(
                inrange_ix, AmGroup, P, E, G, rot_mat, FILES, RunPar, WS
            )

        # Done with Hamiltonian, on to dipoles! (while a separate subroutine
        # would be nice, the Jansen method requires the Estats... don't want
        # to re calculate those... therefore, the dipoles are in here!)
        dipole, dip_r = CalcDipole(
            AmGroup, P, E, G, rot_mat, FILES, RunPar, WS)

        if any(
            outtype in RunPar.output_type for outtype in [
                "Ham", "Dip"
            ]
        ):
            WS.Dipoles[AmGroup, :] = dipole
            WS.TDCgen_r[AmGroup, :] = dip_r

        if "Ram" in RunPar.output_type:
            WS.Raman[AmGroup, :] = CalcRaman(
                AmGroup, P, E, G, rot_mat, FILES, RunPar, WS)


def CalcEstat(inrange_ix, AmGroup, FILES, RunPar, WS):
    """
    Computes the electric properties on this amide group. Takes all residues
    within range, and considers their influence (coulomb only).

    Depending on the choice of map, a different set of atoms is considered
    'local'. Also, the exact method to calculate the influence depends.
    If the map requires the electric gradients, a different (more expensive)
    function is called opposed to when they're not needed.
    C speedups are available for the routines called here, and used when
    applicable.
    """
    if WS.OscID[AmGroup] == 0 and RunPar.map_choice == "Tokmakoff":
        local_ix = WS.ix.local_tok[AmGroup]
    else:
        local_ix = WS.ix.local_dict[AmGroup]

    if RunPar.use_c_lib:
        inrange_ix_c = AIM_SO.difference_c(inrange_ix, local_ix, FILES.clib)
    else:
        inrange_ix = AIM_SO.difference_KvA(inrange_ix, local_ix)

    # 0=C, 1=O, 3=N, 4=H
    P = np.zeros((6), dtype='float32')
    E = np.zeros((6, 3), dtype='float32')
    G = np.zeros((6, 6), dtype='float32')

    # datatype management
    if RunPar.use_c_lib:
        Natoms = np.int32(len(inrange_ix_c))
    else:
        inrange_ix_np = np.array(inrange_ix).astype('int32')
        inrange_pos, inrange_char = AIM_DC.posfinder(
            inrange_ix_np, WS.positions, WS.charges)

    if not WS.OscID[AmGroup] == 0:  # it's not BB amide!
        ex_map = RunPar.ExtraMaps[WS.OscID[AmGroup]]
        calcG = ex_map.use_G
        relevant_j = ex_map.relevant_j
    else:
        if WS.PrePro[AmGroup] == 1:
            calcGforEstat = WS.AllMaps.Emap.Pro.useG
            calcGforDip = WS.AllMaps.Dmap.Pro.useG
            relevant_j = WS.AllMaps.Emap.Pro.relevant_j
        else:
            calcGforEstat = WS.AllMaps.Emap.Gen.useG
            calcGforDip = WS.AllMaps.Dmap.Gen.useG
            relevant_j = WS.AllMaps.Emap.Gen.relevant_j

        if not calcGforDip and not calcGforEstat:
            calcG = False
        else:
            calcG = True

    if calcG:
        if RunPar.use_c_lib:
            for _, j in enumerate(relevant_j):
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                ref_pos_c = np.ctypeslib.as_ctypes(ref_pos)
                P, E, G = CalcFieldGrad_c(ref_pos_c, inrange_ix_c, Natoms, P,
                                          E, G, j, FILES.clib, WS)
        else:
            for _, j in enumerate(relevant_j):
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                P, E, G = CD_CalcFieldGrad(ref_pos, inrange_pos, inrange_char,
                                           P, E, G, j, WS.halfbox,
                                           WS.boxdims)

    else:
        if RunPar.use_c_lib:
            for _, j in enumerate(relevant_j):
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                ref_pos_c = np.ctypeslib.as_ctypes(ref_pos)
                P, E = CalcField_c(ref_pos_c, inrange_ix_c, Natoms, P, E, j,
                                   FILES.clib, WS)
        else:
            for _, j in enumerate(relevant_j):
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                P, E = CD_CalcField(ref_pos, inrange_pos, inrange_char, P, E,
                                    j, WS.halfbox, WS.boxdims)

    return P, E, G


def UpdateEstat(inrange_ix, AmGroup, P, E, G, rot_mat, FILES, RunPar, WS):
    """
    Most notable: Tokmakoff frequencies, but Jansen dipoles. Because they use
    a different set of 'locals', the field is slightly different. Instead of
    completely recalculating the field, it just adds/subtracts the influence
    of the changed locals.
    Also, some atoms weren't required for tokmakoff, but are for jansen: the
    properties of those are also now calculated so they're present for the
    dipole calculation!
    """
    # no updates are required when using Torii dipoles!
    if RunPar.Dipole_choice == "Torii":
        return P, E, G

    # From here: Jansen dipoles!

    if WS.PrePro[AmGroup] == 1:
        relevant_j = WS.AllMaps.Emap.Pro.relevant_j
        calcG = WS.AllMaps.Dmap.Pro.useG
        relevant_j_dipole = WS.AllMaps.Dmap.Pro.relevant_j
    else:
        relevant_j = WS.AllMaps.Emap.Gen.relevant_j
        calcG = WS.AllMaps.Dmap.Gen.useG
        relevant_j_dipole = WS.AllMaps.Dmap.Gen.relevant_j

    relevant_j_new = [j for j in relevant_j_dipole if j not in relevant_j]

    # if Tokmakoff map was used, the data for the calculated atoms has to be
    # corrected for the change in local_ix
    if RunPar.map_choice == "Tokmakoff":
        Pn = np.zeros((6), dtype='float32')
        En = np.zeros((6, 3), dtype='float32')
        Gn = np.zeros((6, 6), dtype='float32')

        local_ix_tok = WS.ix.local_tok[AmGroup]
        local_ix = WS.ix.local_dict[AmGroup]

        local_added = [ix for ix in local_ix[:] if ix not in local_ix_tok[:]]
        local_remov = [ix for ix in local_ix_tok[:] if ix not in local_ix[:]]

        if len(local_added) != 0:
            local_added_np = np.array(local_added).astype('int32')
            inrange_pos, inrange_char = AIM_DC.posfinder(
                local_added_np, WS.positions, WS.charges)

            # added means that they are now local. Their contributions ARE now
            # in PEG, but need to be removed from there.
            inrange_char *= -1

            for j in relevant_j:
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                if np.sum(G[j, :]) == 0:
                    Pn, En = CD_CalcField(
                        ref_pos, inrange_pos, inrange_char, Pn, En, j,
                        WS.halfbox, WS.boxdims)
                else:
                    Pn, En, Gn = CD_CalcFieldGrad(
                        ref_pos, inrange_pos, inrange_char, Pn, En, Gn, j,
                        WS.halfbox, WS.boxdims)

        if len(local_remov) != 0:
            local_remov_np = np.array(local_remov).astype('int32')
            inrange_pos, inrange_char = AIM_DC.posfinder(
                local_remov_np, WS.positions, WS.charges)

            for j in relevant_j:
                ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                if np.sum(G[j, :]) == 0:
                    Pn, En = CD_CalcField(
                        ref_pos, inrange_pos, inrange_char, Pn, En, j,
                        WS.halfbox, WS.boxdims)
                else:
                    Pn, En, Gn = CD_CalcFieldGrad(
                        ref_pos, inrange_pos, inrange_char, Pn, En, Gn, j,
                        WS.halfbox, WS.boxdims)

    else:
        if len(relevant_j_new) == 0:
            return P, E, G
        else:
            Pn = np.zeros((6), dtype='float32')
            En = np.zeros((6, 3), dtype='float32')
            Gn = np.zeros((6, 6), dtype='float32')

    # if the dipole map requires other atoms than the Emap, calculate them also
    if len(relevant_j_new) != 0:
        local_ix = WS.ix.local_dict[AmGroup]
        if RunPar.use_c_lib:
            inrange_ix_c = AIM_SO.difference_c(
                inrange_ix, local_ix, FILES.clib)
        else:
            inrange_ix = AIM_SO.difference_KvA(inrange_ix, local_ix)

        if RunPar.use_c_lib:
            Natoms = np.int32(len(inrange_ix_c))
        else:
            inrange_ix_np = np.array(inrange_ix).astype('int32')
            inrange_pos, inrange_char = AIM_DC.posfinder(
                inrange_ix_np, WS.positions, WS.charges)

        if calcG:
            if RunPar.use_c_lib:
                for j in relevant_j_new:
                    ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                    ref_pos_c = np.ctypeslib.as_ctypes(ref_pos)
                    Pn, En, Gn = CalcFieldGrad_c(
                        ref_pos_c, inrange_ix_c, Natoms, Pn, En, Gn, j,
                        FILES.clib, WS)
            else:
                for j in relevant_j_new:
                    ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                    Pn, En, Gn = CD_CalcFieldGrad(
                        ref_pos, inrange_pos, inrange_char, Pn, En, Gn, j,
                        WS.halfbox, WS.boxdims)
        else:
            if RunPar.use_c_lib:
                for j in relevant_j_new:
                    ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                    ref_pos_c = np.ctypeslib.as_ctypes(ref_pos)
                    Pn, En = CalcField_c(
                        ref_pos_c, inrange_ix_c, Natoms, Pn, En, j, FILES.clib,
                        WS)
            else:
                for j in relevant_j_new:
                    ref_pos = WS.positions[WS.AllOscGroups[AmGroup, 6+j], :]
                    Pn, En = CD_CalcField(
                        ref_pos, inrange_pos, inrange_char, Pn, En, j,
                        WS.halfbox, WS.boxdims)

    # now, both changes of existing atoms and the addition of new atoms has
    # been done. Finalize the changes by merging with the original P, E, G.
    # (of course, after rotating to the correct orientation)
    En, Gn = AIM_MF.EG_rotation(En, Gn, rot_mat, np.arange(6))

    Pn *= RunPar.bohr2ang
    En *= RunPar.bohr2ang2
    Gn *= RunPar.bohr2ang3

    P += Pn
    E += En
    G += Gn

    return P, E, G


@njit
def CD_CalcField(ref_pos, atpos, atchar, P, E, j, halfbox, boxdims):
    """
    The numba way of calculating the electric potential and field on a given
    atom due to the influence of all residues in range.
    """
    for atom in range(atpos.shape[0]):
        diff = AIM_MF.PBC_diff(ref_pos, atpos[atom, :], halfbox, boxdims)
        dist = AIM_MF.vec3_len(diff)
        prefac = atchar[atom] / (dist*dist*dist)
        P[j] += atchar[atom] / dist
        for i in range(3):
            E[j, i] += diff[i] * prefac

    return P, E


def CalcField_c(ref_pos_c, inrange_ix_c, Natoms, P, E, j, clib, WS):
    """
    The c way of calculating the electric potential and field on a given atom
    due to the influence of all residues in range.
    """
    clib.CalcField(ref_pos_c, inrange_ix_c, WS.positions_c, WS.charges_c,
                   WS.halfbox_c, WS.boxdims_c, Natoms, WS.out_c)
    P[j] += WS.out_c[0]
    E[j, :] += WS.out_c[1:4]
    return P, E


@njit
def CD_CalcFieldGrad(ref_pos, atpos, atchar, P, E, G, j,
                     halfbox, boxdims):
    """
    The numba way of calculating the electric potential, field and gradient on
    a given atom due to the influence of all residues in range.
    """
    for atom in range(atpos.shape[0]):
        diff = AIM_MF.PBC_diff(ref_pos, atpos[atom, :], halfbox, boxdims)
        idist2 = 1/AIM_MF.dotprod(diff, diff)
        idist = np.sqrt(idist2)
        prefac = atchar[atom] * (idist*idist2)
        prefac2 = 3.0 * prefac * idist2

        diffX = diff[0]
        diffY = diff[1]
        diffZ = diff[2]

        P[j] += atchar[atom] * idist

        E[j, 0] += diffX * prefac
        E[j, 1] += diffY * prefac
        E[j, 2] += diffZ * prefac

        G[j, 0] += prefac - (diffX * diffX * prefac2)
        G[j, 1] += prefac - (diffY * diffY * prefac2)
        G[j, 2] += prefac - (diffZ * diffZ * prefac2)
        G[j, 3] -= diffX * diffY * prefac2
        G[j, 4] -= diffX * diffZ * prefac2
        G[j, 5] -= diffY * diffZ * prefac2

    return P, E, G


def CalcFieldGrad_c(ref_pos_c, inrange_ix_c, Natoms, P, E, G, j,
                    clib, WS):
    """
    The C way of calculating the electric potential, field and gradient on a
    given atom due to the influence of all residues in range.
    """
    clib.CalcFieldGrad(ref_pos_c, inrange_ix_c, WS.positions_c, WS.charges_c,
                       WS.halfbox_c, WS.boxdims_c, Natoms, WS.out_c)
    P[j] += WS.out_c[0]
    E[j, :] += WS.out_c[1:4]
    G[j, :] += WS.out_c[4:]
    return P, E, G


def CalcHamDiagRaw(OscGroup, P, E, G, FILES, RunPar, WS):
    """
    Takes the computed electric properties and applies a map on them to
    calculate the expected frequency for this group.
    """
    if WS.OscID[OscGroup] == 0:
        freq, rot_mat = calc_freq_BB(
            OscGroup, P, E, G, RunPar, WS)
    else:
        ex_map = RunPar.ExtraMaps[WS.OscID[OscGroup]]
        freq, rot_mat = (
            ex_map.functions["calc_freq"](
                P, E, G, OscGroup, FILES, RunPar, WS, ex_map)
        )
    if "Ham" in RunPar.output_type:
        WS.Hamiltonian[OscGroup, OscGroup] = freq

    return rot_mat


def calc_freq_BB(AmGroup, P, E, G, RunPar, WS):
    """
    Calculating the frequency of the group using just the electric properties
    for amide groups in the protein backbones.
    """
    # the calculated fields have length units in angstrom, but the maps are
    # designed for length units in bohr. /ang2bohr = *bohr2ang
    P *= RunPar.bohr2ang
    E *= RunPar.bohr2ang2
    G *= RunPar.bohr2ang3

    # now, it is time to orient the field properties such that they align
    # with the map direction

    Am_C = WS.AllOscGroups[AmGroup, 6]
    Am_O = WS.AllOscGroups[AmGroup, 7]
    Am_N = WS.AllOscGroups[AmGroup, 9]

    # determine the coordinate system with respect to the amide bond
    Xvec = AIM_MF.PBC_diff(WS.positions[Am_O, :], WS.positions[Am_C, :],
                           WS.halfbox, WS.boxdims)
    Xvec /= AIM_MF.vec3_len(Xvec)
    Yvec = AIM_MF.project(Xvec, AIM_MF.PBC_diff(WS.positions[Am_N, :],
                          WS.positions[Am_C, :], WS.halfbox, WS.boxdims))
    Yvec /= AIM_MF.vec3_len(Yvec)
    Zvec = AIM_MF.crossprod(Xvec, Yvec)

    # maybe, in theory, not necessary, but not that expensive in grand scheme,
    # plus, improves accuracy!
    Zvec /= AIM_MF.vec3_len(Zvec)

    rot_mat = np.empty((3, 3), dtype='float32')
    rot_mat[0, :] = Xvec
    rot_mat[1, :] = Yvec
    rot_mat[2, :] = Zvec

    if WS.PrePro[AmGroup] == 1:
        omega = WS.AllMaps.Emap.Pro.omega
        Pmap = WS.AllMaps.Emap.Pro.Pmap
        Emap = WS.AllMaps.Emap.Pro.Emap
        Gmap = WS.AllMaps.Emap.Pro.Gmap
        relevant_j = WS.AllMaps.Emap.Pro.relevant_j_arr
    else:
        omega = WS.AllMaps.Emap.Gen.omega
        Pmap = WS.AllMaps.Emap.Gen.Pmap
        Emap = WS.AllMaps.Emap.Gen.Emap
        Gmap = WS.AllMaps.Emap.Gen.Gmap
        relevant_j = WS.AllMaps.Emap.Gen.relevant_j_arr

    E, G = AIM_MF.EG_rotation(E, G, rot_mat, relevant_j)

    freq = AIM_MF.apply_map(P, E, G, omega, Pmap, Emap, Gmap, relevant_j)

    return freq, rot_mat


def nnmap(AmGroup, function_warning, logfilename, RunPar, WS):
    """
    In case of the protein backbone, the calculated frequencies after applying
    a map might not be accurate enough. This applies a nearest-neighbour
    correction to the frequency based on the ramachandran angles with the
    neighbours.
    """
    NNarray = WS.AllOscGroups[AmGroup, :]

    # Calculate influence of N terminal neighbour (earlier in the chain)
    if WS.NtermN[AmGroup]:
        # figure out which of the seven maps to use
        mapname = cis_or_trans(NNarray, 1, RunPar, WS)
        NNmap = getattr(WS.AllMaps.NN, "NtermShift" + mapname)

        delta, function_warning = calcdelta(
            NNarray, NNmap, 0, function_warning, logfilename, RunPar, WS)
        WS.Hamiltonian[AmGroup, AmGroup] += delta

    # Calculate influence of C terminal neighbour (further in the chain)
    if WS.CtermN[AmGroup]:
        # figure out which of the seven maps to use
        mapname = cis_or_trans(NNarray, 2, RunPar, WS)
        NNmap = getattr(WS.AllMaps.NN, "CtermShift" + mapname)

        delta, function_warning = calcdelta(
            NNarray, NNmap, 6, function_warning, logfilename, RunPar, WS)
        WS.Hamiltonian[AmGroup, AmGroup] += delta

    return function_warning


def cis_or_trans(NNarray, i, RunPar, WS):
    """
    Given two backbone amide oscillators, decides what kind of NN map to use.
    This depends on whether either of these is or isn't a pre-pro. When there
    is a proline, the N becomes chiral, and this matters.
    """
    # no prolines! easy!
    if (WS.resnames[NNarray[i*6]] != "PRO"
            and WS.resnames[NNarray[i*6+3]] != "PRO"):
        maptype = ""

    # for now, the pro-pro case is treated as if it is pro-gly. in the future,
    # a pro-pro map should be made!
    else:
        r1 = AIM_MF.PBC_diff(WS.positions[NNarray[(i-1)*6+1], :],
                             WS.positions[NNarray[(i-1)*6], :],
                             WS.halfbox, WS.boxdims)  # CO vec
        r2 = AIM_MF.PBC_diff(WS.positions[NNarray[(i-1)*6+4], :],
                             WS.positions[NNarray[(i-1)*6+3], :],
                             WS.halfbox, WS.boxdims)  # NH vec

        # determine cis-trans, based on whether vectors are more aligned, or
        # more opposed.

        # (both pro-pro(should be treated as pro-gly) and pro-gly)
        if WS.resnames[NNarray[i*6]] == "PRO":
            # In the original code, Pro-Pro is actually treated as Gly-Pro
            if (WS.resnames[NNarray[i*6+3]] == "PRO"
                    and RunPar.replicate_orig_AIM):
                BondType = "GP"
            else:
                BondType = "PG"
        else:
            BondType = "GP"

        if BondType == "PG":
            resnum = WS.resnums[NNarray[i*6]]
            LD = WS.residues.LD[resnum]
            if LD < 0:
                Dpro = True
            else:
                Dpro = False

            if AIM_MF.dotprod(r1, r2) < 0:
                if Dpro:
                    maptype = "_transDPro_transGly"
                else:
                    maptype = "_transPro_transGly"
            else:
                if Dpro:
                    maptype = "_cisDPro_transGly"
                else:
                    maptype = "_cisPro_transGly"
        else:
            if AIM_MF.dotprod(r1, r2) < 0:
                maptype = "_transGly_transPro"
            else:
                maptype = "_cisGly_transPro"

    return maptype


def calcdelta(NNarray, NNmap, i, function_warning, logfilename, RunPar, WS):
    """
    Takes the two amide groups and calculates the ramachandran angles between
    them. Then, it plugs these angles in to the supplied map (chosen by the
    cis_or_trans function), obtaining the value delta which is returned.
    """
    dim = 13
    space = 30

    # calculate the ramachandran angles
    phi_ang = AIM_MF.dihedral(
        WS.positions[NNarray[i], :], WS.positions[NNarray[i+3], :],
        WS.positions[NNarray[i+5], :], WS.positions[NNarray[i+6], :],
        WS.halfbox, WS.boxdims) * 180 / np.pi
    psi_ang = AIM_MF.dihedral(
        WS.positions[NNarray[i+3], :], WS.positions[NNarray[i+5], :],
        WS.positions[NNarray[i+6], :], WS.positions[NNarray[i+9], :],
        WS.halfbox, WS.boxdims) * 180 / np.pi

    # binning: see between which values (only lower bound is found as upper
    # is just 1 further) the angles fall
    phi_N = int((phi_ang + 180)//space)
    psi_N = int((psi_ang + 180)//space)
    if phi_N == dim - 1:
        phi_N = dim-2
    if psi_N == dim - 1:
        psi_N = dim-2

    if phi_N >= 0 and phi_N < dim-1 and psi_N >= 0 and psi_N < dim-1:
        # determine lower and higher bound
        x1l = phi_N * space - 180
        x2l = psi_N * space - 180

        y1 = NNmap[psi_N, phi_N]
        y2 = NNmap[psi_N+1, phi_N]
        y3 = NNmap[psi_N+1, phi_N+1]
        y4 = NNmap[psi_N, phi_N+1]

        u = (phi_ang - x1l)/space
        t = (psi_ang - x2l)/space

        # bilinear interpolation!
        delta = (1-u)*(1-t)*y1 + (1-u)*t*y2 + u*t*y3 + u*(1-t)*y4

    else:
        ErrorText = (
            "Ill defined ramachandran angles between residue "
            + str(WS.resnums[NNarray[i]])
            + " and residue " + str(WS.resnums[NNarray[i+3]])
            + " in frame" + str(WS.framenum)
            + ". The nearest neighbour shift will be set to zero.")
        function_warning = AIM_PC.warning(ErrorText, False, function_warning,
                                          logfilename, RunPar)
    return delta, function_warning


def CalcDipole(OscGroup, P, E, G, rot_mat, FILES, RunPar, WS):
    """
    Calculates the dipole moment of this oscillator, along with the position
    of this dipole moment (the latter is required for calculating couplings
    with this dipole moment)
    """
    # CAUTION!!! PEG are still in the format they're left in by
    # CalcHamDiagRaw!
    if WS.OscID[OscGroup] == 0:
        mi, ri = CalcDipole_BB(OscGroup, P, E, G, rot_mat,
                               RunPar, WS)
    else:
        ex_map = RunPar.ExtraMaps[WS.OscID[OscGroup]]
        mi, ri = ex_map.functions["calc_dipole"](
            P, E, G, OscGroup, rot_mat, FILES, RunPar, WS, ex_map
        )
    return mi, ri


def CalcDipole_BB(AmGroup, P, E, G, rot_mat, RunPar, WS):
    """
    Calculates the dipole moment and position of it for amide groups in the
    protein backbone.
    """
    Am_C = WS.AllOscGroups[AmGroup, 6]
    Am_O = WS.AllOscGroups[AmGroup, 7]

    COvec = AIM_MF.PBC_diff(WS.positions[Am_O, :], WS.positions[Am_C, :],
                            WS.halfbox, WS.boxdims)
    COvec /= AIM_MF.vec3_len(COvec)

    if RunPar.Dipole_choice == "Torii":
        Am_N = WS.AllOscGroups[AmGroup, 9]

        CNvec = AIM_MF.PBC_diff(WS.positions[Am_N, :], WS.positions[Am_C, :],
                                WS.halfbox, WS.boxdims)
        CNvec /= AIM_MF.vec3_len(CNvec)

        if RunPar.replicate_orig_AIM:
            magnitude = 2.73
        else:
            magnitude = 0.276

        mi = Dipole_Torii(COvec, CNvec, magnitude)
        ri = WS.positions[Am_C, :] + 0.665*COvec + 0.258*CNvec

    else:
        mi_loc = np.zeros((3))
        mi = np.zeros((3))

        if WS.PrePro[AmGroup] == 1:
            typename = "Pro"
        else:
            typename = "Gen"

        Emap = getattr(WS.AllMaps.Dmap, typename)
        for xyz in range(3):
            mi_loc[xyz] = AIM_MF.apply_map(
                P, E, G, Emap.omega[xyz],
                Emap.Pmap[xyz], Emap.Emap[xyz],
                Emap.Gmap[xyz],
                Emap.relevant_j_arr
            )

        mi = np.dot(mi_loc, rot_mat)

        ri = WS.positions[Am_C, :] + 0.868*COvec

    return mi, ri


@njit
def Dipole_Torii(COvec, CNvec, magnitude):
    """
    Calculates the dipole moment using the Torii method.
    """
    # itheta = 1/0.17632698  ## 1/tan(10 degrees expressed in radians)
    itheta = 5.6712818196

    dri = 0.665*COvec + 0.258*CNvec

    dridri = AIM_MF.dotprod(dri, dri)
    COvecdri = AIM_MF.dotprod(COvec, dri)
    mi = dri - (COvecdri + np.sqrt(dridri - COvecdri*COvecdri)*itheta)*COvec

    mi /= AIM_MF.vec3_len(mi)
    mi *= magnitude

    return mi


def CalcRaman(OscGroup, P, E, G, rot_mat, FILES, RunPar, WS):
    """
    Calculates the Raman tensor of this oscillator.
    Raman tensor in form of a vector with 6 values (xx xy xz yy yz zz)
    """
    if WS.OscID[OscGroup] == 0:
        ra = CalcRaman_BB(OscGroup, P, E, G, rot_mat, WS)
    else:
        ex_map = RunPar.ExtraMaps[WS.OscID[OscGroup]]
        ra = ex_map.functions["calc_raman"](
            OscGroup, P, E, G, rot_mat, FILES, RunPar, WS, ex_map)

    return ra


def CalcRaman_BB(OscGroup, P, E, G, rot_mat, WS):
    """
    Calculates the Raman tensor for amide groups in the protein backbone.
    """
    # get the basis from the rotation matrix
    COvec = rot_mat[0, :]
    CNvec = rot_mat[1, :]
    Zvec = rot_mat[2, :]

    theta = 0.5934119  # 34*np.pi/180

    RXvec = np.cos(theta) * COvec - np.sin(theta) * CNvec
    RYvec = np.sin(theta) * COvec + np.cos(theta) * CNvec

    Rvec = (
        AIM_MF.tensorprod(RXvec)*20 +
        AIM_MF.tensorprod(RYvec)*4 +
        AIM_MF.tensorprod(Zvec)
    )

    return Rvec


def FindAtomPos(WS):
    """
    For each oscillator, takes the position of a key atom. These positions are
    printed to a position file which can be used for CD calculations.
    """
    for grInd, atInd in enumerate(WS.ix.AllAtomPos):
        WS.AtomPos[grInd, :] = WS.positions[atInd, :]
