# standard lib imports
import os
import subprocess
import sys

# 3rd party lib imports
import MDAnalysis as MDA
import numpy as np

# my lib imports
import sourcecode.DataConversion as AIM_DC
import sourcecode.MathFunctions as AIM_MF
import sourcecode.PrintCommands as AIM_PC
import sourcecode.PhysicsFunctions as AIM_PF
import sourcecode.SourceMapreader as AIM_SM
import sourcecode.SetOperations as AIM_SO


class Universe:
    def __init__(self, FILES, RunPar):
        """
        The first basic initialization of the universe/md system. Creates a
        MDA Universe class (saved as Universe.Universe). Checks whether the
        box of the mda system is cubic, whether it has an interger (preferably
        zero) charge, and creates the basic property arrays.
        """
        # generate universe
        self.Universe = gen_universe(FILES, RunPar)

        # test if box has right angles (and save as variable)
        self.rightangled = box_angle_check(FILES, RunPar, self.Universe)

        # generate the arrays storing all the system's properties
        self.SetProperties(FILES, RunPar)

        # test if the whole system has zero charge. Quits for non-integer
        # charge, raises warning (but continues) for integer non-zero charge.
        charge_check(FILES, RunPar, self.charges)

        # get the maps that are actually needed (mainly, get the correct E and
        # D map in place)
        self.AllMaps = AIM_SM.Needed_Maps(RunPar, FILES.AllMaps)

    def SetProperties(self, FILES, RunPar):
        """
        Creates the arrays storing the properties of all the atoms in the
        system. If required, also creates the c-friendly version of these
        arrays.
        """
        self.atnums = self.Universe.atoms.ix
        self.atnames = self.Universe.atoms.names
        self.resnums = self.Universe.atoms.resnums
        self.resnames = self.Universe.atoms.resnames
        self.positions = self.Universe.atoms.positions
        self.masses = self.Universe.atoms.masses
        self.charges = self.Universe.atoms.charges
        self.types = self.Universe.atoms.types
        self.segids = self.Universe.atoms.segids

        self.natoms = np.int32(self.resnums.shape[0])

        # CHARMM (and others?) resets resnum with each residue
        self.AbsResnums(FILES, RunPar)
        self.nres = np.int32(self.resnums[-1] + 1)

        self.boxdims = self.Universe.dimensions[:3].astype('float32')
        self.halfbox = self.boxdims/2
        self.halfbox = self.halfbox.astype('float32')

        self.residues = ResidueFinder(FILES.logfilename, RunPar, self)
        self.ix = IXFinder(self.residues)

        # Gromacs has molnums appointed -> dimer protein = mols,
        # every solvent/lipid residue = molnum.

        if not RunPar.atom_based_chainID:
            try:
                self.molnums = self.Universe.atoms.molnums
            except MDA.exceptions.NoDataError:
                self.FindMolnums(FILES, RunPar)
        else:
            self.FindMolnums_atID(FILES, RunPar)

        self.residues.protein_init(FILES.logfilename, RunPar, self)

        RunPar.SetOPLS(self.types[0][:4] == "opls")

        if RunPar.use_c_lib:
            self.charges = self.charges.astype('float32')
            self.charges_c = np.ctypeslib.as_ctypes(self.charges)
            self.masses = self.masses.astype('float32')
            self.masses_c = np.ctypeslib.as_ctypes(self.masses)
            self.res_COM_c = np.ctypeslib.as_ctypes(
                np.zeros((self.nres * 3), dtype='float32'))

    def FindMolnums(self, FILES, RunPar):
        if RunPar.atom_based_chainID:
            self.FindMolnums_atID(FILES, RunPar)
        else:
            self.FindMolnums_SegID(RunPar)

    def FindMolnums_SegID(self, RunPar):
        """
        CHARMM support. Determines molnums using segment and residue numbers.
        """
        self.molnums = np.zeros(self.atnums.shape, dtype='int32')

        protinseg = {}
        segid = ""
        for atomnum, resname in enumerate(self.resnames):
            prevsegid = segid
            segid = self.segids[atomnum]
            if segid == prevsegid:
                if protinseg[segid]:
                    continue
                if resname in RunPar.resnames_Protein:
                    protinseg[segid] = True
            else:
                if resname in RunPar.resnames_Protein:
                    protinseg[segid] = True
                else:
                    protinseg[segid] = False

        molnum = -1
        prevsegid = ""
        prevresnum = -1
        for atomnum, segid in enumerate(self.segids):
            resnum = self.resnums[atomnum]
            if protinseg[segid]:
                if segid != prevsegid:
                    molnum += 1
            else:
                if resnum != prevresnum:
                    molnum += 1

            self.molnums[atomnum] = molnum
            prevsegid = segid
            prevresnum = self.resnums[atomnum]

    def FindMolnums_atID(self, FILES, RunPar):
        """
        Amber support. Determines molnums using atom names and bonds.
        """
        self.molnums = np.zeros(self.atnums.shape, dtype='int32')

        u = self.Universe

        prev_resnum = -1
        molnum = -1
        stop_mol = False
        resnums_in_current_molnum = set()
        for atnum, resnum in enumerate(self.resnums):

            # If not the first atom of new residue:

            # we already know what to do with this residue!
            if resnum == prev_resnum:
                self.molnums[atnum] = molnum
                continue

            # -------------

            # If it is the first atom of a residue:

            # general preparation.
            resname = self.resnames[atnum]
            prev_resnum = resnum  # prev_resnum won't be used below

            # if it's not a protein, every residue is a molecule
            if resname not in RunPar.resnames_Protein:
                # In case a protein will come after this non-prot residue:
                # That alone should mean that a new molecule starts.
                stop_mol = True
                molnum += 1
                self.molnums[atnum] = molnum
                continue

            # ---------------

            # From here: a new protein residue:

            # if we detected that the previous residue MUST've proven that
            # a new chain is starting:
            # this is done either by not being a protein, or by that protein
            # chain showing it's ended.
            if stop_mol:
                stop_mol = False
                molnum += 1
                self.molnums[atnum] = molnum
                resnums_in_current_molnum = {resnum}
                continue

            # some residue names MUST be the beginning of a new chain,
            # just like some atom names MUST indicate the beginning of a new
            # chain
            if (
                resname in ["FOR"]
                or
                any(
                    self.atnames[ix] in RunPar.atnames["TermH"]
                    for ix in range(
                        self.residues.start[resnum],
                        self.residues.fin[resnum] + 1
                    )
                )
            ):
                molnum += 1
                self.molnums[atnum] = molnum
                resnums_in_current_molnum = {resnum}
                continue

            # Similarly, some residue names MUST be the end of a chain, just
            # like some atom names MUST indicate the beginning of a new chain
            if(
                resname in ["ETA", "GL2"]
                or
                any(
                    self.atnames[ix] in RunPar.atnames["TermO"]
                    for ix in range(
                        self.residues.start[resnum],
                        self.residues.fin[resnum] + 1
                    )
                )
            ):
                stop_mol = True
                self.molnums[atnum] = molnum
                continue

            # ----------

            # This is all we can do using labels (atname, atnum, resname,
            # resnum, etc) to find out whether a chain started/ended.

            # From here:
            # - mid-chain residues
            # - cyclic-chain residues
            # - some other start/end to a chain not identified using TermH and
            #   TermO atoms, or resname.

            # So, only possible identification left: bonds! Find the C of this
            # chain, and retrieve the resnum of the N that it is bonded to.
            # if it's a new resnum for this chain -> this was mid-chain!
            # if it's a resnum that occured earlier in this chain -> must be
            # cyclic, and this is the last residue!

            # firstly: only before this loop, resnums_in_current_molnum is
            # defined as an empty set. In ALL other instances, it is defined
            # containing the resnum that triggered formation. Therefore, it
            # being empty means a new chain!
            if len(resnums_in_current_molnum) == 0:
                molnum += 1
                self.molnums[atnum] = molnum
                resnums_in_current_molnum = {resnum}
                continue

            # now, lets find out if our C is bound to a N!

            first_at = u.atoms[atnum]  # an MDA atomtype.

            # We need to find the C in the same residue as this atom. Lets use
            # indices!
            # Gromacs, Amber -> all resnums are unique
            # Charmm -> recurring resnums, but each has different segid.

            cann_resnum = str(first_at.resnum)
            cann_segid = str(first_at.segid)

            # MDA byres should also work?
            this_residue = u.atoms[
                self.residues.start[resnum]:self.residues.fin[resnum]+1
            ]
            Cgr = this_residue.select_atoms(
                "name " + " or name ".join(RunPar.atnames["C"])
            )

            # C is an atomgroup!
            if len(Cgr) != 1:
                # This should not be possible. 'normal' protein residues (i.e.
                # amino acids) should have exactly 1 generic C (the C=O).
                # If this is a special one (one that terminates the chain),
                # it should be programmed in (like FOR).
                # If there's multiple, there's some other geometry issue...
                ErrorText = (
                    "\nUh-oh! While trying to identify the protein chains in "
                    "the system, a residue was found that was expected to be "
                    "a part of the protein chain, but which doesn't contain "
                    "exactly one carbon atom of the type that binds to the "
                    "nitrogen of the next residue!\nTo be precise, this "
                    "residue "
                )
                if len(Cgr) == 0:
                    ErrorText += (
                        "doesn't contain any of these."
                    )
                else:
                    ErrorText += (
                        "contains " + str(len(Cgr)) + " of these."
                    )
                ErrorText += (
                    " Please check the manual for what to do now! Some info "
                    "on this residue:"
                    "\nresidue number: " + str(resnum) +
                    "\nresidue name: " + str(resname) +
                    "\ncannonical residue number: " + cann_resnum +
                    "\ncannonical segment ID: " + cann_segid +
                    "\ndesired type(s) of atom: " + ", ".join(
                        RunPar.atnames["C"]) +
                    "\n\n Quitting!"
                )
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)

            # the atom that we're expecting to bind to the N of the next
            # residue
            Cat = Cgr[0]
            C_atnum = Cat.ix

            first_resnum = min(resnums_in_current_molnum)
            first_res = u.atoms[
                self.residues.start[first_resnum]:
                self.residues.fin[first_resnum]+1,
            ]
            next_res = u.atoms[
                self.residues.start[resnum+1]:
                self.residues.fin[resnum+1]+1
            ]
            subgroup = first_res + this_residue + next_res

            Ngr = subgroup.select_atoms(
                "( name " + " or name ".join(RunPar.atnames["N"]) +
                " ) and bonded index " + str(C_atnum)
            )

            if len(Ngr) > 1:
                ErrorText = (
                    "\nUh-oh! While trying to identify the protein chains in "
                    "the system, a residue was found that has two N-terminal "
                    "neighbours. Please consult the manual! Some info on this "
                    "residue:"
                    "\nresidue number: " + str(resnum) +
                    "\nresidue name: " + str(resname) +
                    "\ncannonical residue number: " + cann_resnum +
                    "\ncannonical segment ID: " + cann_segid +
                    "\ndesired type(s) of atom: " + ", ".join(
                        RunPar.atnames["N"]) +
                    "\n\n Quitting!"
                )
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)

            elif len(Ngr) == 0:
                # last residue of chain!
                stop_mol = True
                self.molnums[atnum] = molnum
                continue

            # now, there is definitely a neighbour.
            Nat = Ngr[0]
            N_resnum = Nat.resnum

            if N_resnum in resnums_in_current_molnum:
                # 'end' of a cyclic chain!
                stop_mol = True
            else:
                # not yet the end of this chain
                resnums_in_current_molnum.add(N_resnum)
            self.molnums[atnum] = molnum

    def AbsResnums(self, FILES, RunPar):
        """
        CHARMM support. Makes sure that every residue number is unique, and
        that the first one has ID 0.
        """
        prevresnum = -1
        writeresnum = -1
        for atomnum, ix in enumerate(self.atnums):
            if atomnum != ix:
                ErrorText = (
                    "The atom number of the atom at position "
                    + str(atomnum)
                    + " does not match the position in the list! Quitting!")
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)
            resnum = self.resnums[atomnum]
            if resnum != prevresnum:
                writeresnum += 1
            self.resnums[atomnum] = writeresnum
            prevresnum = resnum

    def InitSys(self, FILES, RunPar):
        """
        Initializes the system. Detects the geometry, and creates some lookup
        tables that are heavily used later on. No more big tables are created
        after this function.
        """

        logfilename = FILES.logfilename

        # creates AllAmGroups array, N/CtermN, BBSC, PrePro
        self.osc_finder(logfilename, RunPar)
        self.residues.COM_init(RunPar, self)

        if RunPar.use_c_lib:
            self.residues.c_ify(self)
            self.ix.c_ify(self)
            self.out_c = np.ctypeslib.as_ctypes(
                np.zeros((10), dtype='float32'))

        if "Ham" in RunPar.output_type:
            self.InitCouplingArrays(RunPar)
        else:
            if "Dip" in RunPar.output_type:
                self.TDCgen_r = np.empty((self.res_desired_len, 3))

        self.characteristics = Characterizer(self)
        self.ix.LocalFinder(RunPar, self)

    def osc_finder(self, logfilename, RunPar):
        """
        This is the function responsible for finding all oscillators in the
        system. After this function has run, the AllOscGroups array exists and
        contains the indices of all atoms used in the simulation.
        """
        # FinishCalc needs to know something about the system, but nothing in
        # that regard is saved... everything it needs to know is gathered in
        # this dict.
        self.fincalcprints = {}

        dict_of_BBCOatoms, dict_of_cycles = self.osc_finder_BB(
            logfilename, RunPar)

        for molnum, arr in dict_of_BBCOatoms.items():
            self.fincalcprints["nBBraw" + str(molnum)] = arr.shape[0]

        dict_of_extra_osc = self.osc_finder_extra(RunPar)

        self.fincalcprints["OscID_present"] = [
            k for k in dict_of_extra_osc.keys()]

        self.SystemSize_new(
            dict_of_BBCOatoms, dict_of_cycles, dict_of_extra_osc, logfilename,
            RunPar)

        self.OscArr_builder(
            dict_of_BBCOatoms, dict_of_cycles, dict_of_extra_osc,
            RunPar
        )

    def osc_finder_BB(self, logfilename, RunPar):
        """
        Finds all amide oscillators contained in the protein backbones.
        CAUTION: this function assumes a single protein chain covers precisely
        1 molum index. Reversely, one molecule (molnum index) containing a
        protein chain contains only and all amino acids belonging to that
        chain, and no others.
        Creates an array of shape(n, 6) for each protein chain.
        """
        # find out how many protein chains there are
        protein_molnums = []
        for resnum in self.residues.resnums_protein:
            molnum = self.molnums[self.residues.start[resnum]]
            if molnum not in protein_molnums:
                protein_molnums.append(molnum)
        n_mols = len(protein_molnums)
        self.fincalcprints["molnum_prot_chains"] = protein_molnums
        AIM_PC.vprintl(
            1, ["\nTotal amount of detected chains:", n_mols],
            logfilename, RunPar)

        dict_of_cycles = {}
        dict_of_BBCOatoms = {}

        # This should help solve issues where protein chains are not correctly
        # identified, perhaps because something else has a protein-type
        # residue name. AIM is designed for proteins, but maybe, this way, it
        # could also function for systems without them?
        if "AmideBB" not in RunPar.oscillators:
            return dict_of_BBCOatoms, dict_of_cycles

        for molnum in protein_molnums:
            molres = []

            # find out which residue numbers are in this molecule
            for resnum, atnum in enumerate(self.residues.start):
                if self.molnums[atnum] == molnum:
                    molres.append(resnum)
                elif self.molnums[atnum] > molnum:
                    break

            # Check if indeed, this molecule contains nothing but protein:
            for resnum in molres:
                if (
                    self.resnames[self.residues.start[resnum]]
                    not in RunPar.resnames_Protein
                ):
                    ErrorText = (
                        "\nMolecule number " + str(molnum) + " has been found"
                        " to contain both residues that are expected to be "
                        "part of a protein chain, and residues that are not. "
                        "If your system contains a protein, please refer to "
                        "the manual. If not, don't include 'AmideBB' (without"
                        " quotation marks) in the list of oscillators "
                        "specified in the input file. Quitting!"
                    )
                    AIM_PC.warning(ErrorText, True, 0, logfilename, RunPar)

            # find out whether this protein chain is cyclic
            dict_of_cycles[molnum] = True

            # look at atoms in first residue of molecule:
            for atomnum in range(
                self.residues.start[molres[0]],
                self.residues.fin[molres[0]]+1
            ):
                if self.resnames[atomnum] in ["FOR"]:
                    dict_of_cycles[molnum] = False
                    break
                if self.atnames[atomnum] in RunPar.atnames["TermH"]:
                    dict_of_cycles[molnum] = False
                    break

            # look at atoms in last residue of molecule:
            for atomnum in range(
                self.residues.start[molres[-1]],
                self.residues.fin[molres[-1]]+1
            ):
                if self.resnames[atomnum] in ["ETA", "GL2"]:
                    dict_of_cycles[molnum] = False
                    break
                if self.atnames[atomnum] in RunPar.atnames["TermO"]:
                    dict_of_cycles[molnum] = False
                    break

            # print the findings on whether is protein chain is cyclic
            n_res = len(molres)
            if dict_of_cycles[molnum]:
                AIM_PC.vprintl(1, ["Chain", molnum,
                                   "is identified as a cyclic chain of", n_res,
                                   "residues"], logfilename, RunPar)
            else:
                AIM_PC.vprintl(1, ["Chain", molnum,
                                   "is identified as a linear chain of", n_res,
                                   "residues"], logfilename, RunPar)

            # array is 1 shorter than the amount of residues, unless the chain
            # is cyclic
            relevant_atoms = np.empty(
                (n_res-1+dict_of_cycles[molnum], 6),
                dtype='int32'
            )
            found_atoms = np.ones(
                (n_res - 1 + dict_of_cycles[molnum], 6),
                dtype='int32'
            )

            # now, fill the array with atoms:
            for locresnum, resnum in enumerate(molres):
                resname = self.residues.resnames[resnum]

                for atnum in range(
                    self.residues.start[resnum],
                    self.residues.fin[resnum]+1
                ):
                    atname = self.atnames[atnum]

                    # if it is not the last residue of linear chain
                    if locresnum != n_res-1 or dict_of_cycles[molnum]:
                        # save the atoms. 0->C_x, 1->O_x, 2->CA_x
                        if atname in RunPar.atnames["C"]:
                            relevant_atoms[locresnum, 0] = atnum
                            found_atoms[locresnum, 0] = 0
                        elif atname in RunPar.atnames["O"]:
                            relevant_atoms[locresnum, 1] = atnum
                            found_atoms[locresnum, 1] = 0
                        elif atname in RunPar.atnames["CA"]:
                            relevant_atoms[locresnum, 2] = atnum
                            found_atoms[locresnum, 2] = 0
                        elif resname == "FOR":
                            if atname == "H":
                                relevant_atoms[locresnum, 2] = atnum
                                found_atoms[locresnum, 2] = 0

                    # if it is not the first residue of linear chain
                    if locresnum != 0 or dict_of_cycles[molnum]:
                        # save the atoms. 3->N_x+1, 4->H_x+1 (or CD_x+1),
                        # 5->CA_x+1
                        if atname in RunPar.atnames["N"]:
                            relevant_atoms[locresnum-1, 3] = atnum
                            found_atoms[locresnum-1, 3] = 0
                        elif atname in RunPar.atnames["H"]:
                            relevant_atoms[locresnum-1, 4] = atnum
                            found_atoms[locresnum-1, 4] = 0
                        elif atname in RunPar.atnames["CA"]:
                            relevant_atoms[locresnum-1, 5] = atnum
                            found_atoms[locresnum-1, 5] = 0
                        elif resname == "ETA":
                            if atname == "C2":
                                relevant_atoms[locresnum-1, 5] = atnum
                                found_atoms[locresnum-1, 5] = 0
                        elif resname == "GL2":
                            if atname == "H1":
                                relevant_atoms[locresnum-1, 4] = atnum
                                found_atoms[locresnum-1, 4] = 0
                            if atname == "H2":
                                relevant_atoms[locresnum-1, 5] = atnum
                                found_atoms[locresnum-1, 5] = 0
                        if resname == "PRO":
                            if atname == "CD":
                                relevant_atoms[locresnum-1, 4] = atnum
                                found_atoms[locresnum-1, 4] = 0
            dict_of_BBCOatoms[molnum] = relevant_atoms

            if np.sum(found_atoms) != 0:
                problem_atoms = np.sum(found_atoms, axis=0)
                problem_dict = {
                    0: "C", 1: "O", 2: "CA", 3: "N", 4: "H", 5: "CA"
                }
                problem_txt = ", ".join([
                    problem_dict[ix] for ix, summed in enumerate(problem_atoms)
                    if summed != 0
                ])
                ErrorText = (
                    "Chain " + str(molnum) + " containing residues " +
                    str(molres[0]) + " - " + str(molres[-1]) +
                    " contains " + str(np.sum(found_atoms)) +
                    " atoms which are needed for amide groups, but could " +
                    "not be identified. These atoms are of type " + problem_txt
                    + ". If you are missing many atoms of a single type, "
                    "please check if your atoms have the same name as reported"
                    " in the used atnames file. Quitting!"
                )
                AIM_PC.warning(ErrorText, True, 0, logfilename, RunPar)

        return dict_of_BBCOatoms, dict_of_cycles

    def osc_finder_extra(self, RunPar):
        """
        Finds all oscillators encoded in any extramaps file. Creates an array
        of shape (n, 6) for each type of oscillator.
        """
        raw_relatoms_dict = {}

        # We want to create a similar dict to that of osc_finder_BB. Here, the
        # keys will be the ID of the groups, the value the array containing
        # all oscillators. Every oscillator gets a row.

        allmaps = []
        for extramap in RunPar.ExtraMaps.values():
            if extramap.name in RunPar.oscillators:
                allmaps.append(extramap)

        for resnum, resname in enumerate(self.residues.resnames):
            for extramap in allmaps:
                if resname not in extramap.makeup:
                    continue
                osc_atoms = self.find_osc_atoms(
                    extramap.makeup[resname], resnum)
                if len(osc_atoms) == 0:
                    continue
                if extramap.ID in raw_relatoms_dict:
                    raw_relatoms_dict[extramap.ID].extend(osc_atoms)
                else:
                    raw_relatoms_dict[extramap.ID] = osc_atoms

        relevant_atoms_dict = {}
        for ID, relatlist in raw_relatoms_dict.items():
            n_osc = len(relatlist)
            defarr = np.empty((n_osc, 6), dtype='int32')
            for ix, subarr in enumerate(relatlist):
                defarr[ix, :] = subarr

            if RunPar.ExtraMaps[ID].name in RunPar.oscillators:
                relevant_atoms_dict[ID] = defarr

        return relevant_atoms_dict

    def find_osc_atoms(self, allowed_configs, resnum):
        """
        If a given residue has the name that a map is looking for, test whether
        the correct atoms are present also. Return all sets of 6 atoms that
        fit for this residue/map combination.
        """
        osc_atoms = []
        # there might be multiple ways to define this map for this residue
        for allowed_config in allowed_configs:
            check = 0
            # go check all atoms!
            atoms = np.empty((6), dtype='int32')
            for atnum in range(
                self.residues.start[resnum],
                self.residues.fin[resnum]+1
            ):
                atname = self.atnames[atnum]
                for i in range(6):
                    if atname in allowed_config[i]:
                        atoms[i] = atnum
                        check += 10 ** (5-i)
            if check == 111111:
                osc_atoms.append(atoms)
        return osc_atoms

    def SystemSize_new(
        self, dict_of_BBCOatoms, dict_of_cycles, dict_of_extra_osc,
        logfilename, RunPar
    ):
        """
        Determines how many oscillators should be considered for this system
        in total, for this calculation.
        """
        n_BB_CO_groups = 0
        for molnum, array in dict_of_BBCOatoms.items():
            n_BB_CO_groups += array.shape[0]
            if not dict_of_cycles[molnum]:
                lastresnum = self.resnums[array[-1, 3]]
                # although the residue number was allowed according to the
                # black/whitelists, this residue number does not have an amide
                # bond associated with it -> remove from the desired resnums
                # list
                if lastresnum in self.residues.desired:
                    self.residues.desired.remove(lastresnum)

        n_other_groups = 0
        for array in dict_of_extra_osc.values():
            n_other_groups += array.shape[0]

        self.res_desired_len = np.int32(
            len(self.residues.desired) * ("AmideBB" in RunPar.oscillators)
            + n_other_groups
        )

        if self.res_desired_len == 0:
            ErrorText = (
                "\nSelection criteria for the oscillators are too strict, no "
                "oscillators match the requirements, and therefore, the "
                "frequency of none of them should be calculated. Quitting!"
            )
            AIM_PC.warning(ErrorText, True, 0, logfilename, RunPar)

    def OscArr_builder(
        self, dict_of_BBCOatoms, dict_of_cycles, dict_of_extra_osc, RunPar
    ):
        """
        Creates the most important arrays with information on the oscillators.
        These store the indices of all atoms that make up the considered
        oscillators, whether these oscillators have neighbouring oscillators,
        or a proline neighbour (applicable only for protein backbone amide
        groups), and what the ID number of their type is.
        """
        # create all storage

        # for each osc, one line with all important atom numbers
        self.AllOscGroups = np.empty((self.res_desired_len, 18), dtype='int32')

        # for each osc, whether it has a C-terminal neighbour (BBAm only)
        self.CtermN = np.empty((self.res_desired_len), dtype=bool)

        # for each osc, whether it has a N-terminal neighbour (BBAm only)
        self.NtermN = np.empty((self.res_desired_len), dtype=bool)

        # for each osc, it's ID number
        self.OscID = np.empty((self.res_desired_len), dtype='int32')

        # for each osc, whether its a prePro or not (BBAm only)
        self.PrePro = np.empty((self.res_desired_len), dtype="int32")

        counter = 0
        if "AmideBB" in RunPar.oscillators:
            # loop over all chains in the backbone dict
            for molnum, array in dict_of_BBCOatoms.items():
                Ngroups = array.shape[0]
                used_osc_in_chain = 0
                for AmGroup in range(Ngroups):

                    # for the white/blacklist -> if the residue number is not
                    # selected, skip it!
                    resnum = self.resnums[array[AmGroup, 0]]
                    if resnum not in self.residues.desired:
                        continue
                    # The middle 6 entries are for the amide group itself.
                    self.AllOscGroups[counter, 6:12] = array[AmGroup, :]
                    self.OscID[counter] = 0

                    # prePro detection - no.9 is the N atom of AmGroup
                    self.PrePro[counter] = int(
                        self.resnames[self.AllOscGroups[counter, 9]] == "PRO"
                    )

                    # Neighbour finding
                    if AmGroup == 0:
                        if dict_of_cycles[molnum]:
                            # if first residue of cyclic chain -> N-term
                            # neighbour is the last group of the same chain
                            self.NtermN[counter] = True
                            self.AllOscGroups[counter, 0:6] = array[-1, :]
                        else:
                            # if first residue of linear chain -> no N-term
                            # neighbour present
                            self.NtermN[counter] = False
                            self.AllOscGroups[counter, 0:6] = 999999999
                    else:
                        # if not first residue, there must be a N-term
                        # neighbour (these are BB only)
                        self.NtermN[counter] = True
                        self.AllOscGroups[counter, 0:6] = array[AmGroup-1, :]

                    if AmGroup == Ngroups-1:
                        if dict_of_cycles[molnum]:
                            # if last group of a cyclic chain -> C-term
                            # neighbour is the first group of the same chain
                            self.CtermN[counter] = True
                            self.AllOscGroups[counter, 12:] = array[0, :]
                        else:
                            # if last group of linear chain -> no C-term
                            # neighbour
                            self.CtermN[counter] = False
                            self.AllOscGroups[counter, 12:] = 999999999
                    else:
                        # if not last residue, there must be a C-term neighbour
                        self.CtermN[counter] = True
                        self.AllOscGroups[counter, 12:] = array[AmGroup+1, :]

                    counter += 1
                    used_osc_in_chain += 1
                self.fincalcprints["nBBused" + str(molnum)] = used_osc_in_chain

        for IDnum, array in dict_of_extra_osc.items():
            Ngroups = array.shape[0]
            for OscGroup in range(Ngroups):
                self.AllOscGroups[counter, 6:12] = array[OscGroup, :]
                self.OscID[counter] = IDnum
                self.NtermN[counter] = False
                self.AllOscGroups[counter, 0:6] = 999999999
                self.CtermN[counter] = False
                self.AllOscGroups[counter, 12:] = 999999999
                self.PrePro[counter] = 0
                counter += 1

        if RunPar.use_c_lib:
            self.AllOscGroups_c = AIM_DC.ctype2d(self.AllOscGroups, 'int32')
            self.PrePro_c = np.ctypeslib.as_ctypes(self.PrePro)

    def InitCouplingArrays(self, RunPar):
        """
        Creates the np arrays that will store all kinds of information
        regarding coupling.
        """
        if (RunPar.coupling_choice == "TCC"
                or RunPar.NN_coupling_choice == "TCC"):
            self.TCC_v = np.empty((self.res_desired_len, 18))
            if RunPar.use_c_lib:
                self.TCC_v_c = AIM_DC.ctype2d(self.TCC_v, 'float32')
        if (RunPar.coupling_choice == "TDCKrimm"
                or RunPar.NN_coupling_choice == "TDCKrimm"):
            self.TDCKr_r = np.empty((self.res_desired_len, 3))
            self.TDCKr_m = np.empty((self.res_desired_len, 3))
            if RunPar.use_c_lib:
                self.TDCKr_r_c = AIM_DC.ctype2d(self.TDCKr_r, 'float32')
                self.TDCKr_m_c = AIM_DC.ctype2d(self.TDCKr_m, 'float32')
        # why initialize TDCTasumi for TCC-non neighbours? when running in
        # orig_aim mode, the coupling choice works slightly differently... For
        # some reason, the original code only allowed TDC methods for coupling
        # between two non-BB amide groups. If one didn't choose the Krimm
        # method, Tasumi was always applied, even if the non-neighbour choice
        # was set to TCC.
        # The orig_aim check is not executed here. The tasumi arrays are only
        # filled and used in orig_aim mode, and not always!
        if (RunPar.coupling_choice == "TDCTasumi"
                or RunPar.NN_coupling_choice == "TDCTasumi"
                or RunPar.coupling_choice == "TCC"):
            self.TDCTa_r = np.empty((self.res_desired_len, 3))
            self.TDCTa_m = np.empty((self.res_desired_len, 3))
            if RunPar.use_c_lib:
                self.TDCTa_r_c = AIM_DC.ctype2d(self.TDCTa_r, 'float32')
                self.TDCTa_m_c = AIM_DC.ctype2d(self.TDCTa_m, 'float32')

        # for the generic TD coupling between two arbitrary groups
        # There is no need to initialize a TDCgen_m, as it would just equal
        # WS.Dipoles. That one is used instead.
        self.TDCgen_r = np.empty((self.res_desired_len, 3))

    def RunCalc(self, TIMER, FILES, RunPar):
        """
        The master function for the tough part of the calculation: This is
        where the 'for each frame' loop lives. Makes sure everything is
        calculated, and the right things are printed.
        """
        function_warning = 0
        self.PreCalc(TIMER, FILES, RunPar)

        for frame in self.trj[RunPar.start_frame:]:
            # if frame number is too large, quit.
            self.framenum = frame.frame
            if self.framenum < RunPar.start_frame:
                continue
            if self.framenum == RunPar.start_frame:
                TIMER.StartCalc()
            self.relframenum = self.framenum - RunPar.start_frame

            self.PrintFrameNum(TIMER, FILES.logfilename, RunPar)

            if self.framenum >= RunPar.end_frame:
                break

            if self.CheckFrameNum(TIMER, FILES.logfilename, RunPar):
                break  # manages the frame number prints

            # positions update, spheresize check
            function_warning, breaker = self.UpdateData(
                FILES.logfilename, RunPar, frame, function_warning)
            if breaker:
                break

            self.COM_update(FILES, RunPar)
            self.CreateHam(RunPar)

            # if anything needs to happen to the map on a per-frame basis.
            for exmap in RunPar.ExtraMaps.values():
                exmap.functions["pre_frame"](FILES, RunPar, self, exmap)
            for coupmap in RunPar.CouplingMaps.values():
                coupmap.functions["pre_frame"](FILES, RunPar, self, coupmap)

            AIM_PC.vprint(
                4, ("Frame initialization complete, starting calculating "
                    "Hamiltonian"),
                FILES.logfilename, RunPar)
            AIM_PF.CalcHam(FILES, RunPar, self)
            AIM_PC.vprint(
                4,
                "finished calculating Hamiltonian, writing ham to file",
                FILES.logfilename, RunPar)

            for exmap in RunPar.ExtraMaps.values():
                exmap.functions["post_frame"](FILES, RunPar, self, exmap)
            for coupmap in RunPar.CouplingMaps.values():
                coupmap.functions["post_calc"](FILES, RunPar, self, coupmap)

            FILES.WriteOutput(RunPar, self)

            self.lastCalcFrame = self.framenum

    def PreCalc(self, TIMER, FILES, RunPar):
        """
        Takes care of some business (including prints) that needs to happen
        just before starting the heavy calculations.
        """
        self.NSA_atdict, self.NSA_resdict, self.NSA_lendict = {}, {}, {}
        AIM_PC.vprint(
            2,
            "\n" + FILES.script_version + " Starting iterating over frames",
            FILES.logfilename, RunPar)

        TIMER.AfterInit()

        if "Ham" in RunPar.output_type:
            self.GenCoupTypeArr(FILES, RunPar)

        # if people want to do any pre-calc work with their own maps
        for exmap in RunPar.ExtraMaps.values():
            exmap.functions["pre_calc"](FILES, RunPar, self, exmap)
        for coupmap in RunPar.CouplingMaps.values():
            coupmap.functions["pre_calc"](FILES, RunPar, self, coupmap)

        self.trj = self.Universe.trajectory

        if RunPar.end_frame >= len(self.trj):
            RunPar.end_frame = len(self.trj)
        self.max_loop = 0

    def GenCoupTypeArr(self, FILES, RunPar):
        """
        Creates the array that stores what coupling map should be used for any
        given pair of oscillators.
        """
        self.coup_type_array = np.zeros((
            self.res_desired_len, self.res_desired_len))
        self.coup_used = {}

        def_coup_dict = {
            "None": 0,
            "TDCKrimm": 2,
            "TDCTasumi": 3,
            "TCC": 4,
            "Tasumi": 101,
            "GLDP": 102
        }
        NN_ID = def_coup_dict[RunPar.NN_coupling_choice]
        non_NN_ID = def_coup_dict[RunPar.coupling_choice]

        all_pairs = RunPar.Coupling_all_pairs
        for AmGroupi in range(self.res_desired_len):
            CAin = self.AllOscGroups[AmGroupi, 2]
            CAic = self.AllOscGroups[AmGroupi, 14]
            OscIDi = self.OscID[AmGroupi]
            for AmGroupj in range(AmGroupi):
                CAj = self.AllOscGroups[AmGroupj, 8]
                OscIDj = self.OscID[AmGroupj]
                pair = (min(OscIDi, OscIDj), max(OscIDi, OscIDj))
                groups = [AmGroupi, AmGroupj]

                # if the pair has been defined especially using a .cmap file
                if pair in all_pairs:
                    ID = all_pairs[pair]
                    self.coup_type_array[AmGroupi, AmGroupj] = ID
                    self.coup_type_array[AmGroupj, AmGroupi] = ID
                    if ID not in self.coup_used:
                        self.coup_used[ID] = set()
                    self.coup_used[ID].update(groups)
                    continue

                # if the pair is an amide-only one
                if all(k in [0, 1] for k in pair):
                    # if it's neighbours
                    if CAj == CAin or CAj == CAic:
                        # original version treated opposite ends of cyclic
                        # chain differently: always GLDP
                        if (
                            RunPar.replicate_orig_AIM
                            and
                            CAj == CAic
                        ):
                            tempID = def_coup_dict["GLDP"]
                        else:
                            tempID = NN_ID
                        self.coup_type_array[AmGroupi, AmGroupj] = tempID
                        self.coup_type_array[AmGroupj, AmGroupi] = tempID
                        if tempID not in self.coup_used:
                            self.coup_used[tempID] = set()
                        self.coup_used[tempID].update(groups)
                    # if it's not neighbours
                    else:
                        tempID = non_NN_ID
                        # original version forced TDC method on SC-SC inter-
                        # actions. Pair contains the ID's of the oscillators,
                        # 1 -> AmideSC.
                        if (
                            RunPar.replicate_orig_AIM
                            and
                            all(k == 1 for k in pair)
                            and
                            RunPar.coupling_choice == "TCC"
                        ):
                            tempID = def_coup_dict["TDCTasumi"]
                        self.coup_type_array[AmGroupi, AmGroupj] = tempID
                        self.coup_type_array[AmGroupj, AmGroupi] = tempID
                        if tempID not in self.coup_used:
                            self.coup_used[tempID] = set()
                        self.coup_used[tempID].update(groups)
                    continue

                # Now, the only remaining groups are not amide-only (so AIM
                # doesn't know what to do with them), and they have not been
                # specified by the user using .cmap files either. Now, only
                # two options; apply default dipole-dipole coupling, or
                # don't do anything (don't calculate coupling at all). Which
                # depends on the apply_dd_coupling parameter.

                if OscIDi == OscIDj:  # they're the same kind of group
                    if RunPar.apply_dd_coupling in ["All", "Same"]:
                        self.coup_type_array[AmGroupi, AmGroupj] = 1
                        self.coup_type_array[AmGroupj, AmGroupi] = 1
                        if 1 not in self.coup_used:
                            self.coup_used[1] = set()
                        self.coup_used[1].update(groups)
                else:  # they're two different kinds of groups
                    if RunPar.apply_dd_coupling in ["All", "Different"]:
                        self.coup_type_array[AmGroupi, AmGroupj] = 1
                        self.coup_type_array[AmGroupj, AmGroupi] = 1
                        if 1 not in self.coup_used:
                            self.coup_used[1] = set()
                        self.coup_used[1].update(groups)
        for coupID, group_set in self.coup_used.items():
            self.coup_used[coupID] = np.array(
                sorted(list(group_set)), dtype='int32')
        if RunPar.use_c_lib:
            self.coup_used_c = {
                k: np.ctypeslib.as_ctypes(v)
                for k, v in self.coup_used.items()}

    def PrintFrameNum(self, TIMER, logfilename, RunPar):
        """
        This function prints the current frame number and ETA.
        """
        max_loop_count = 0
        relframenum_copy = self.relframenum
        while True:
            if relframenum_copy < 10:
                if relframenum_copy == 1:
                    printlevel = 1
                else:
                    printlevel = 2

                TIMER.Update()
                frames_todo = RunPar.end_frame - self.framenum

                if self.framenum >= RunPar.end_frame:
                    pretext = "Stopping on frame "
                else:
                    pretext = "Starting on frame "

                try:
                    time_avgperframe = ((TIMER.update - TIMER.startheavycalc)
                                        / self.relframenum)
                except ZeroDivisionError:
                    AIM_PC.vprint(
                        printlevel, pretext + str(self.framenum),
                        logfilename, RunPar)
                else:
                    time_remaining = int(time_avgperframe * frames_todo)
                    AIM_PC.vprint(
                        printlevel, pretext + str(self.framenum) + "."
                        + (10-max_loop_count) * " " + "Elapsed time = "
                        + AIM_DC.timestring(int(TIMER.update))
                        + "   Estimated remaining time = "
                        + AIM_DC.timestring(time_remaining),
                        logfilename, RunPar)
                if max_loop_count > self.max_loop:
                    self.max_loop = max_loop_count
                    AIM_PC.vprint(2, " ", logfilename, RunPar)
                break
            elif relframenum_copy % 10 == 0:
                relframenum_copy /= 10
                max_loop_count += 1
            else:
                if RunPar.Verbose == 4 or RunPar.Verbose_log == 4:
                    TIMER.Update()
                    frames_todo = RunPar.end_frame - self.framenum

                    if self.framenum >= RunPar.end_frame:
                        pretext = "Stopping on frame "
                    else:
                        pretext = "Starting on frame "

                    try:
                        time_avgperframe = ((TIMER.update
                                             - TIMER.startheavycalc)
                                            / self.relframenum)
                    except ZeroDivisionError:
                        AIM_PC.vprint(
                            4, pretext + str(self.framenum),
                            logfilename, RunPar)
                    else:
                        time_remaining = int(time_avgperframe * frames_todo)
                        AIM_PC.vprint(
                            4, pretext + str(self.framenum) + "."
                            + (10-max_loop_count) * " " + "Elapsed time = "
                            + AIM_DC.timestring(int(TIMER.update))
                            + "   Estimated remaining time = "
                            + AIM_DC.timestring(time_remaining),
                            logfilename, RunPar)
                break

    def CheckFrameNum(self, TIMER, logfilename, RunPar):
        """
        Every 25th frame, it is determined whether the calculation should
        continue, or has to be aborted (because of time restrictions)
        """
        returnval = False
        if self.relframenum % 25 == 0:
            TIMER.Update()
            try:
                time_avgperframe = ((TIMER.update - TIMER.startheavycalc)
                                    / self.relframenum)
            except ZeroDivisionError:
                pass
            else:
                time_left = RunPar.max_time - TIMER.update
                time_perbatch = time_avgperframe * 25
                if time_left < 2 * time_perbatch:
                    AIM_PC.vprint(
                        0, "\nNot enough time left to start next batch. "
                        + "Finalizing calculation. To continue next time, "
                        + "use 'start_frame " + str(self.framenum) + "'.",
                        logfilename, RunPar)
                    returnval = True

        return returnval

    def UpdateData(self, logfilename, RunPar, frame, function_warning):
        """
        Updates the system info for this new frame. The positions array has to
        be refreshed, and it must be checked whether the new box still fits the
        requested SphereSize.
        """
        # build the new position arrays (these have to be updated every frame)
        self.positions = self.Universe.atoms.positions
        self.boxdims = frame.dimensions[:3].astype('float32')
        self.halfbox = self.boxdims/2
        self.halfbox = self.halfbox.astype('float32')
        if RunPar.use_c_lib:
            self.positions_c = AIM_DC.ctype2d(self.positions, 'float32')
            self.boxdims_c = np.ctypeslib.as_ctypes(self.boxdims)
            self.halfbox_c = np.ctypeslib.as_ctypes(self.halfbox)

        returnval = False
        # checking SphereSize:
        if self.framenum == RunPar.start_frame:
            if RunPar.SphereSize > min(self.halfbox):
                ErrorText = ("Selected SphereSize is too large, it may not be "
                             "larger than half of the smallest box dimension! "
                             "Quitting!")
                AIM_PC.warning(ErrorText, True, 0, logfilename, RunPar)
        else:
            if RunPar.SphereSize > 1.1 * min(self.halfbox):
                ErrorText = ("Selected SphereSize is too large, it should "
                             "never become larger than half of the smallest "
                             "box dimension! As now, it is more than 10% "
                             "larger than that... Stopping the run at frame "
                             "number " + str(self.framenum) + "!")
                function_warning_new = AIM_PC.warning(
                    ErrorText, False, function_warning, logfilename, RunPar)
                returnval = True
            elif RunPar.SphereSize > min(self.halfbox):
                ErrorText = ("Selected SphereSize is too large, it should "
                             "never become larger than half of the smallest "
                             "box dimension! While initially correct, the box "
                             "has shrunk enough for it to become too small for"
                             " the selected spheresize. Affected frame: "
                             + str(self.framenum))
                function_warning_new = AIM_PC.warning(
                    ErrorText, False, function_warning, logfilename, RunPar)
        if "function_warning_new" not in locals():
            function_warning_new = function_warning

        return function_warning_new, returnval

    def COM_update(self, FILES, RunPar):
        """
        As the positions changed, so did the COM. This function recalculates
        the COM.
        """
        if RunPar.NSA_toggle:
            # every NSA_nframes'th frame, do this:
            if (self.framenum % RunPar.NSA_nframes == 0
                    or self.relframenum == 0):
                # determine Centres of Mass of all residues we're interested
                # in (part of the type selection).
                # the resulting Res_COM is a dictionary of all COMs, where
                # Key = residue.resnum, Value = np array of size 3, containing
                # the x, y and z coordinates of the centre of mass.

                self.Res_COM = AIM_PF.calc_COM("Sel", FILES, RunPar, self)

                all_inrange = np.zeros((self.nres), dtype='int32')
                # for each AmGroup that we consider:
                for AmGroup in range(self.res_desired_len):
                    resnum = self.resnums[self.AllOscGroups[AmGroup, 6]]
                    NSA_sphere_resnums = AIM_PF.AGsorterNSA(
                        resnum, FILES, RunPar, self)

                    # All inrange sees whether any given residue is considered
                    # for any AmGroup. If after the for AmGroup range a value
                    # is still zero, it is too far away from all calculated
                    # groups.
                    for resindex in NSA_sphere_resnums:
                        all_inrange[resindex] = 1

                    # save the shortlist in a dictionary.
                    if RunPar.use_c_lib:
                        self.NSA_lendict[resnum] = np.int32(len(
                            NSA_sphere_resnums))
                        self.NSA_resdict[resnum] = np.ctypeslib.as_ctypes(
                            NSA_sphere_resnums)
                    else:
                        self.NSA_resdict[resnum] = NSA_sphere_resnums

                # create a list which contains all atom/residue numbers that
                # have a '1' in the all_inrange array.
                self.all_inrange_nres = np.sum(all_inrange)
                self.all_inrange_res = np.empty(
                    (self.all_inrange_nres), dtype='int32')
                counter = 0
                for i in range(self.nres):
                    if all_inrange[i] == 1:
                        self.all_inrange_res[counter] = i
                        counter += 1

                # data type management:
                if RunPar.use_c_lib:
                    self.all_inrange_res_c = np.ctypeslib.as_ctypes(
                        self.all_inrange_res)

            else:
                # this is why the all-residue lists were made -> don't have to
                # calculate COM of a residue if it is not in range anyways.
                self.Res_COM = AIM_PF.calc_COM("NSA", FILES, RunPar, self)

        # NSA_toggle == False! In that case, just calculate COM, and continue.
        else:
            self.Res_COM = AIM_PF.calc_COM("Sel", FILES, RunPar, self)

    def CreateHam(self, RunPar):
        """
        Initialize the Hamiltonian, Dipole and AtomPos arrays.
        """

        if "Ham" in RunPar.output_type:
            self.Hamiltonian = np.zeros(
                (self.res_desired_len, self.res_desired_len), dtype='float32')

        # the dipoles are also used to calculate couplings!
        if any(i in RunPar.output_type for i in ["Ham", "Dip"]):
            self.Dipoles = np.zeros((self.res_desired_len, 3), dtype='float32')
        if "Pos" in RunPar.output_type:
            self.AtomPos = np.zeros((self.res_desired_len, 3), dtype='float32')

        if "Ram" in RunPar.output_type:  # Raman#
            self.Raman = np.zeros((self.res_desired_len, 6), dtype='float32')

    def FinishCalc(self, TIMER, FILES, RunPar):
        """
        Takes care of all the prints that are done after the calculation
        finishes. Also takes care of finalizing the profiler: its results are
        stored in a file, and if requested, also converted to a png figure.
        """
        AIM_PC.vprint(
            1, "\n\n====================\nCalculation "
            "Finished!\n====================\n",
            FILES.logfilename, RunPar)

        totaloftype = AIM_PC.finprint_composition(FILES, RunPar, self)

        AIM_PC.vprintl(1, ["\n", "\nNotes from maps:"], FILES.logfilename, RunPar)

        for exmap in RunPar.ExtraMaps.values():
            exmap.functions["post_calc"](FILES, RunPar, self, exmap)
        for coupmap in RunPar.CouplingMaps.values():
            coupmap.functions["post_calc"](FILES, RunPar, self, coupmap)

        AIM_PC.finprint_frames(FILES, RunPar, self)

        AIM_PC.finprint_time(FILES, RunPar, self, TIMER)

        AIM_PC.finprint_references(FILES, RunPar, self, totaloftype)

        # if the profiler ran, stop it, and export the collected data
        if RunPar.profiler:
            FILES.pr.create_stats()
            FILES.pr.dump_stats(FILES.proffilename)

        # if so desired, make a pretty graph of the runtime
        if RunPar.pngout:
            strcommand = [
                sys.executable,
                os.path.join(FILES.script_dir, "gprof2dot.py"),
                "-f", "pstats", FILES.proffilename, "-o",
                os.path.join(FILES.script_dir, "temp.out")]
            subprocess.run(strcommand)

            strcommand = [
                "dot", "-Tpng", "-o", FILES.pngfilename,
                os.path.join(FILES.script_dir, "temp.out")]
            subprocess.run(strcommand)

            os.remove(os.path.join(FILES.script_dir, "temp.out"))

        # Aaaaand.... we're done!
        AIM_PC.vprint(
            1, "\n" + "*"*72 +
            "\nIf you read this, the program encountered no problems, and "
            "ran normally!\n" + "*"*72, FILES.logfilename, RunPar)


class ResidueFinder:
    """
    The class that takes care of the residue identification. It saves the
    lowest and highest atnum present in this residue, and, for protein
    backbone residues, whether they are D or L. Finally, any it checks
    whether all residue names are known/expected/familiar. __init__ takes care
    of all of this.
    """
    def __init__(self, logfilename, RunPar, WS):

        self.build_srf(WS, RunPar, logfilename)

    def protein_init(self, logfilename, RunPar, WS):
        """
        Does protein-related identification. For each residue, it is determined
        whether it is L or D. Also, for each residue, it checks whether the
        residue should be considered based on the black-/whitelists.
        """
        self.LD = np.zeros((WS.nres))
        self.desired = []

        if len(self.resnums_protein) == 0:
            ProtRes_first = -1
            ProtRes_last = -1
        else:
            ProtRes_first = self.resnums_protein[0]
            ProtRes_last = self.resnums_protein[-1]
        for resnum, resname in enumerate(self.resnames):
            if resnum < ProtRes_first:
                continue
            if resnum > ProtRes_last:
                break

            if resnum in self.resnums_protein:
                self.LD[resnum] = self.DLcheck(resnum, logfilename, RunPar, WS)

                if RunPar.Use_AmGroup_selection_criteria:
                    if (RunPar.Use_resname_blacklist
                            and resname in RunPar.resname_blacklist):
                        continue
                    if (RunPar.Use_resname_whitelist
                            and resname not in RunPar.resname_whitelist):
                        continue
                    if (RunPar.Use_resnum_blacklist
                            and resnum in RunPar.resnum_blacklist):
                        continue
                    if (RunPar.Use_resnum_whitelist
                            and resnum not in RunPar.resnum_whitelist):
                        continue

                self.desired.append(resnum)

    def build_srf(self, WS, RunPar, logfilename):
        """
        Create the start, fin, and resname arrays. These arrays have dimension
        nresiudes, not natoms, to make retrieving the residue properties
        considerably easier during calculation.
        """
        self.start = []
        self.resnames = []
        self.fin = []
        self.resnums_influencers = []
        self.resnums_protein = []
        unknowns = []

        last_resnum = -1
        current_resnum = -1

        # loop over all atoms in system
        for atnum, resname in enumerate(WS.resnames):
            last_resnum = current_resnum
            current_resnum = WS.resnums[atnum]

            # detect if this atom is the first of a new residue. If not, we're
            # not interested right this moment, and want to consider the next
            # atom

            if last_resnum == current_resnum:
                pass
            else:
                if last_resnum != -1:
                    self.fin.append(atnum - 1)
                self.resnames.append(resname)
                self.start.append(atnum)

                if resname in RunPar.resnames_Protein:
                    self.resnums_protein.append(current_resnum)

                if resname in RunPar.resnames_influencers:
                    self.resnums_influencers.append(current_resnum)
                else:
                    if (
                        resname not in RunPar.resnames_All
                        and
                        resname not in unknowns
                    ):
                        unknowns.append(resname)

            if atnum == WS.natoms - 1:
                self.fin.append(atnum)

        self.start = np.array(self.start)
        self.fin = np.array(self.fin)

        # if a residue is unknown, it might cause problems later down the line.
        # Therefore, raise error now, and quit!
        self.report_unknowns(unknowns, logfilename, RunPar)

    def report_unknowns(self, unknowns, logfilename, RunPar):
        """
        If unknown residues were found (name not in resnames file), report
        which residue names were unfamiliar.
        """
        if len(unknowns) != 0:
            ErrorText = ("There are " + str(len(unknowns))
                         + " unknown residue types. The first few are:")
            AIM_PC.warning(ErrorText, False, 0, logfilename, RunPar)
            for index, value in enumerate(unknowns):
                if index >= 5:
                    break
                AIM_PC.vprintl(0, [
                    "Residue no.", self.resnames.index(value), "Resname =",
                    value
                ], logfilename, RunPar)
            quit()

    def DLcheck(self, resnum, logfilename, RunPar, WS):
        """
        Determines whether a protein residue is in D- or L-configuration.
        """
        if self.resnames[resnum] in ["GLY", "FOR", "ETA", "GL2"]:
            returnval = 0
        else:
            resstart = self.start[resnum]
            resfin = self.fin[resnum]
            for atom in range(resstart, resfin+1):
                atname = WS.atnames[atom]
                if atname in RunPar.atnames["CA"]:
                    Am_CA = WS.positions[atom, :]
                elif atname in RunPar.atnames["C"]:
                    Am_C = WS.positions[atom, :]
                elif atname in RunPar.atnames["N"]:
                    Am_N = WS.positions[atom, :]
                elif atname == "CB":
                    Am_CB = WS.positions[atom, :]

            for atname in ["CA", "C", "N", "CB"]:
                if "Am_" + atname not in locals():
                    ErrorText = (
                        "Residue number " + str(resnum) + " of type '"
                        + self.resnames[resnum] + "' has no atom of type '"
                        + atname + "'. This type is expected for protein "
                        "residues. Cannot check wether the residue has a D- "
                        "or L- stereocenter, but this will give issues later "
                        "in the calculation, too. \n\nIs this residue really "
                        "supposed to be a part of the protein, or was a "
                        "mistake made when defining the protein using the "
                        "resnames file? If it is supposed to be there, the "
                        "residue is not a standard type, and the code should "
                        "be edited to properly use this residue. Quitting!")
                    AIM_PC.warning(ErrorText, True, 0, logfilename, RunPar)

            Am_C = AIM_MF.PBC_diff(Am_C, Am_CA, WS.halfbox, WS.boxdims)
            Am_N = AIM_MF.PBC_diff(Am_N, Am_CA, WS.halfbox, WS.boxdims)
            Am_CB = AIM_MF.PBC_diff(Am_CB, Am_CA, WS.halfbox, WS.boxdims)

            CxN = AIM_MF.crossprod(Am_C, Am_N)

            if AIM_MF.dotprod(CxN, Am_CB) > 0:
                returnval = -1
            else:
                returnval = 1
        return returnval

    def COM_init(self, RunPar, WS):
        """
        Determines of which groups the COM needs to be calculated (especially
        relevant if the influencers don't just contain all residues)
        """
        nGroups = WS.res_desired_len
        Osc_resnums = set()

        for OscGroupi in range(nGroups):
            resnum = WS.resnums[WS.AllOscGroups[OscGroupi, 6]]
            Osc_resnums.add(resnum)

        Osc_resnums = list(Osc_resnums)
        Osc_resnums.sort()

        COM = AIM_SO.union_KvA(self.resnums_influencers, Osc_resnums)

        self.resnums_protein = np.array(self.resnums_protein)
        self.resnums_influencers = np.array(self.resnums_influencers)
        self.resnums_COM = np.array(COM)

    def c_ify(self, WS):
        """
        Convert data in a c-friendly format where needed.
        """
        self.resnums_COM = self.resnums_COM.astype('int32')
        self.resnums_COM_c = np.ctypeslib.as_ctypes(self.resnums_COM)
        self.resnums_COM_len = np.int32(self.resnums_COM.shape[0])

        self.resnums_influencers = self.resnums_influencers.astype('int32')
        self.resnums_influencers_c = np.ctypeslib.as_ctypes(
            self.resnums_influencers)
        self.resnums_influencers_len = np.int32(
            self.resnums_influencers.shape[0])

        self.start = self.start.astype('int32')
        self.start_c = np.ctypeslib.as_ctypes(self.start)

        self.fin = self.fin.astype('int32')
        self.fin_c = np.ctypeslib.as_ctypes(self.fin)

        self.inrange_c = np.ctypeslib.as_ctypes(
            np.zeros((WS.nres), dtype='int32'))


class IXFinder:
    """
    Deals with some atom-index related issues. Creates/stores the list
    containing all atoms that are part of the protein, as well as the
    dictionary of which atoms are local for each of the oscillators.
    """
    def __init__(self, RESIDUES):
        self.protein = []
        for resnum in RESIDUES.resnums_protein:
            start = RESIDUES.start[resnum]
            fin = RESIDUES.fin[resnum]
            self.protein.extend(range(start, fin+1))

    def c_ify(self, WS):
        self.inrange_long_c = np.ctypeslib.as_ctypes(
            np.zeros((WS.natoms), dtype='int32'))

    def LocalFinder(self, RunPar, WS):
        """
        Creates a dict. For each oscillator (key), what atoms are local
        (value)? A local atom is one that should not be considered when
        calculating the electrostatic potential, field, and gradient.
        """
        self.AllAtomPos = []
        self.local_dict = {}
        self.local_tok = {}

        for AmGroup in range(WS.res_desired_len):
            Am_C = WS.AllOscGroups[AmGroup, 6]
            Am_N = WS.AllOscGroups[AmGroup, 9]
            if RunPar.AtomPos_choice == "C":
                self.AllAtomPos.append(Am_C)
            elif RunPar.AtomPos_choice == "N":
                self.AllAtomPos.append(Am_N)
            elif RunPar.AtomPos_choice == "O":
                self.AllAtomPos.append(WS.AllOscGroups[AmGroup, 7])
            elif RunPar.AtomPos_choice == "D":
                self.AllAtomPos.append(WS.AllOscGroups[AmGroup, 10])

            resnum = WS.resnums[Am_C]
            resnext = WS.resnums[Am_N]

            if WS.OscID[AmGroup] == 0:  # if BB
                local_ix = AIM_SO.sixsort_KvA(WS.AllOscGroups[AmGroup, 6:12])
                local_ix_tok = AIM_SO.sixsort_KvA(
                    WS.AllOscGroups[AmGroup, 6:12])
                if WS.resnames[Am_N] == "PRO":
                    if RunPar.map_choice == "Tokmakoff":
                        # for the Tokmakoff map, only 5 of the current amide
                        # group are locals
                        local_ix_tok = AIM_SO.threesort_KvA(
                            WS.AllOscGroups[AmGroup, 6:9])
                        temp = AIM_SO.twosort_KvA(
                            [WS.AllOscGroups[AmGroup, 9],
                             WS.AllOscGroups[AmGroup, 11]])
                        local_ix_tok = AIM_SO.union_KvA(local_ix_tok, temp)

                    # for other maps, if it is a prePro, don't include the
                    # hydrogen on the CD atom, as well as all 6 of the
                    # current amide group.
                    temp_resnums = WS.characteristics.Pro_HD_resnums
                    temp_atnums = WS.characteristics.Pro_HD_atnums
                    temp_natoms = WS.characteristics.Pro_HD_natoms
                    templist = []
                    for index in range(temp_natoms):
                        if temp_resnums[index] == resnext:
                            templist.append(temp_atnums[index])
                    local_ix = AIM_SO.union_KvA(templist, local_ix)
                sortcritNN = []  # list of residue numbers to consider

                if RunPar.TreatNN:
                    # if it isn't the last of a chain -> Cterm (high index)
                    # influence
                    if WS.CtermN[AmGroup]:
                        # residue number of the group that donated the N to the
                        # C-term neighbour AmGroup
                        Cnn_resnum = WS.resnums[WS.AllOscGroups[AmGroup, 15]]
                        sortcritNN.append(Cnn_resnum)
                        # if a C-term neighbour, add it's 6 atoms to the local
                        # list. However, in the case of ETA residue, DO not add
                        # the atom saved in the last slot (the C2 atom)
                        if WS.atnames[WS.AllOscGroups[AmGroup, -1]] == "C2":
                            temp = AIM_SO.threesort_KvA(
                                WS.AllOscGroups[AmGroup, 12:15])
                            temp.extend(AIM_SO.twosort_KvA(
                                WS.AllOscGroups[AmGroup, 15:17]))
                            local_ix = AIM_SO.union_KvA(temp, local_ix)
                        else:
                            local_ix = AIM_SO.union_KvA(AIM_SO.sixsort_KvA(
                                WS.AllOscGroups[AmGroup, 12:18]), local_ix)

                    temp_resnums = WS.characteristics.Pro_resnums
                    temp_atnums = WS.characteristics.Pro_atnums
                    temp_atnames = WS.characteristics.Pro_atnames
                    temp_natoms = WS.characteristics.Pro_natoms
                    templist = []

                    # if C_neighbor: ( (resnum in [resnum, resnext,
                    # sortcritNN] ) and resname PRO ) and not (resnum
                    # Cnn_resnum and ( name C or name O ) )
                    # else: ( (resnum in [resnum, resnext, sortcritNN] ) and
                    # resname PRO )

                    COnamelist = RunPar.atnames["C"]
                    COnamelist.extend(RunPar.atnames["O"])

                    for index in range(temp_natoms):
                        loop_resnum = temp_resnums[index]
                        if (loop_resnum == resnum
                                or loop_resnum == resnext
                                or loop_resnum in sortcritNN):
                            append = True
                            if (WS.CtermN[AmGroup]
                                    and loop_resnum == Cnn_resnum
                                    and temp_atnames[index] in COnamelist):
                                append = False
                            if append:
                                templist.append(temp_atnums[index])
                    local_ix = AIM_SO.union_KvA(templist, local_ix)

                    # if it isn't the first of a chain -> Nterm (low index)
                    # influence
                    if WS.NtermN[AmGroup]:
                        # for the opls forcefield, the non-CA atoms of the N
                        # term neighbour should be added to locals. Otherwise,
                        # all atoms of the N term neighbour should be added to
                        # locals
                        if RunPar.opls:
                            temp = AIM_SO.union_KvA(
                                AIM_SO.twosort_KvA(
                                    WS.AllOscGroups[AmGroup, :2]),
                                AIM_SO.twosort_KvA(
                                    WS.AllOscGroups[AmGroup, 3:5])
                            )
                            local_ix = AIM_SO.union_KvA(temp, local_ix)
                        else:
                            sortcritNN.append(
                                WS.resnums[WS.AllOscGroups[AmGroup, 0]])
                            local_ix = AIM_SO.union_KvA(AIM_SO.sixsort_KvA(
                                WS.AllOscGroups[AmGroup, :6]), local_ix)

                # if not Tokmakoff map, the hydrogen atoms on the residue's
                # and its neighbours' CA atoms should also be added to the
                # locals!
                temp_resnums = WS.characteristics.Alpha_H_resnums
                temp_atnums = WS.characteristics.Alpha_H_atnums
                temp_natoms = WS.characteristics.Alpha_H_natoms
                templist = []
                for index in range(temp_natoms):
                    loop_resnum = temp_resnums[index]
                    if (loop_resnum == resnum
                            or loop_resnum == resnext
                            or loop_resnum in sortcritNN):
                        templist.append(temp_atnums[index])
                local_ix = AIM_SO.union_KvA(templist, local_ix)

            else:
                # if not BB -> locals can only be within the group itself,
                # nothing else!
                local_ix = (
                    RunPar.ExtraMaps[WS.OscID[AmGroup]].functions[
                        "local_finder"
                    ](
                        WS.AllOscGroups[AmGroup, 6:12]
                    )
                )

            # data-type management! If we're using the C-library, we need this
            # 'locals' list to be in a c-friendly format. Otherwise, just leave
            # it as the list it currently is!
            if RunPar.use_c_lib:
                self.local_dict[AmGroup] = AIM_DC.ctypelist(local_ix, 'int32')
                if RunPar.map_choice == "Tokmakoff":
                    self.local_tok[AmGroup] = AIM_DC.ctypelist(
                        local_ix_tok, 'int32')
            else:
                self.local_dict[AmGroup] = local_ix
                if RunPar.map_choice == "Tokmakoff":
                    self.local_tok[AmGroup] = local_ix_tok


class Characterizer:
    """
    There are some atoms that require special attention. To make searching for
    these atoms easier, they are saved by category in dictionaries stored in
    this class.
    """
    def __init__(self, WS):
        Alpha_H_nums, Pro_nums, Pro_HD_nums = [], [], []

        # loop over all atoms in the protein. If their (residue)name fits one
        # of the categories, their number is appended to the corresponding list
        for index, atnum in enumerate(WS.ix.protein):
            atname = WS.atnames[atnum]
            if atname in ["HA", "HA1", "HA2", "HA3"]:
                Alpha_H_nums.append(atnum)
            elif atname in ["HD1", "HD2"]:
                Pro_HD_nums.append(atnum)
            if WS.resnames[atnum] == "PRO":
                Pro_nums.append(atnum)

        # for each of the lists, build a numpy array with the data, and fill
        # it.
        self.Alpha_H_natoms = len(Alpha_H_nums)
        self.Alpha_H_resnums = np.empty((self.Alpha_H_natoms), dtype='int32')
        self.Alpha_H_atnums = np.empty((self.Alpha_H_natoms), dtype='int32')
        for index, atnum in enumerate(Alpha_H_nums):
            self.Alpha_H_resnums[index] = WS.resnums[atnum]
            self.Alpha_H_atnums[index] = atnum

        self.Pro_natoms = len(Pro_nums)
        self.Pro_resnums = np.empty((self.Pro_natoms), dtype='int32')
        self.Pro_atnums = np.empty((self.Pro_natoms), dtype='int32')
        self.Pro_atnames = np.empty((self.Pro_natoms), dtype=np.dtype("U10"))
        for index, atnum in enumerate(Pro_nums):
            self.Pro_resnums[index] = WS.resnums[atnum]
            self.Pro_atnums[index] = atnum
            self.Pro_atnames[index] = WS.atnames[atnum]

        self.Pro_HD_natoms = len(Pro_HD_nums)
        self.Pro_HD_resnums = np.empty((self.Pro_HD_natoms), dtype='int32')
        self.Pro_HD_atnums = np.empty((self.Pro_HD_natoms), dtype='int32')
        for index, atnum in enumerate(Pro_HD_nums):
            self.Pro_HD_resnums[index] = WS.resnums[atnum]
            self.Pro_HD_atnums[index] = atnum


def gen_universe(FILES, RunPar):
    """
    Tries to generate the universe object from MDA. Catches errors upon
    failure.
    """

    try:
        the_universe = MDA.Universe(FILES.topfile, FILES.trjfile)
    except FileNotFoundError:
        ErrorText = (
            "\nCould not open the topology or trajectory file. Please "
            "make sure the names are correct, and try again. Quitting!")
        AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)
    except TypeError:
        ErrorText = (
            "\nThe given topology/trajectory files are of the wrong type. "
            "Please remember that only some combinations are "
            "currently supported by the program. To see which, check the "
            "manual!\n"
        )
        AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)

    return the_universe


def box_angle_check(FILES, RunPar, the_universe):
    """
    Checks whether the simulation box has right angles. Raises an error if it
    doesn't.
    """
    angs = the_universe.dimensions[3:]
    if angs[0] == 90 and angs[1] == 90 and angs[2] == 90:
        AIM_PC.vprint(3, "\nAll box vectors are perpendicular!",
                      FILES.logfilename, RunPar)
        rightangled = True
    else:
        ErrorText = (
            "\nSubmitted trajectory does not have a cubic, tetragonal or "
            "orthorhombic symmetry. Please make sure all angles of the "
            "simulation box are 90 degrees, as these are the only type of "
            "system AmideImaps can handle. Quitting!")
        AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)
        rightangled = False

    return rightangled


def charge_check(FILES, RunPar, charges_array):
    """
    Checks the total charge of the system. Raises a warning if the charge is
    non-zero, and an error if it is not integer.
    """
    charge_whole_system = charges_array.sum()
    if abs(charge_whole_system) > 0.0001:
        if abs(round(charge_whole_system)-charge_whole_system) > 0.0001:
            ErrorText = ("The total charge of the submitted system is not an "
                         "integer number. Quitting!")
            AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, RunPar)
        else:
            ErrorText = ("The total charge of the system does not equal 0. Be "
                         "aware that this might yield incorrect results!")
            AIM_PC.warning(ErrorText, False, 0, FILES.logfilename, RunPar)
