# standard lib imports
import os

# 3rd party lib imports
import numpy as np

# my lib imports
import sourcecode.DataConversion as AIM_DC
import sourcecode.FileHandler as AIM_FH
import sourcecode.PrintCommands as AIM_PC
import sourcecode.ReferenceManager as AIM_RM


class AllSourceMaps:
    def __init__(self, FILES, RunPar):

        function_warning = 0

        allfuncts = self.get_allfuncts()

        self.available_Emaps = []
        self.available_Dmaps = []

        # find all maps available in the sourcedir folder
        for filename in os.listdir(FILES.sourcedir):

            # skip if the object is not a file
            if not os.path.isfile(os.path.join(FILES.sourcedir, filename)):
                continue

            # skip if the file is not a map file
            if not filename.startswith("Map"):
                continue

            # skip if the file is not of the correct format
            if not filename.endswith(".dat"):
                continue

            # the name without 'Map-' and '.dat'
            # say we have a name 'Map-E_Jansen_backup.dat.
            maptypename = filename[4:-4]  # this becomes 'E_Jansen_backup'
            temp = maptypename.split("_")  # turns it into a '_' separated list
            maptype = temp[0]  # becomes 'E'
            mapname = "_".join(temp[1:])  # becomes 'Jansen_backup'

            full_fname = os.path.join(FILES.sourcedir, filename)
            map = allfuncts[maptype](
                full_fname, maptypename, FILES, RunPar, function_warning)
            function_warning = map.function_warning
            if map.success:
                setattr(self, maptypename, map)
                if maptype == "D":
                    self.available_Dmaps.append(mapname)
                elif maptype == "E":
                    self.available_Emaps.append(mapname)

        exhead = " of class AllSourceMaps"

        AIM_PC.finwarning(
            function_warning, FILES.logfilename, RunPar, extra_header=exhead
        )

    def get_allfuncts(self):
        allfuncts = {
            "C": det_Cmap,
            "D": DE_Map,
            "E": DE_Map
        }
        return allfuncts


class DE_Map:
    def __init__(
        self, full_fname, maptypename, FILES, RunPar, function_warning
    ):
        self.success = True
        self.full_fname = full_fname

        self.maptypename = maptypename
        temp = maptypename.split("_")
        self.maptype = temp[0]
        self.mapname = "_".join(temp[1:])

        self.function_warning = function_warning
        self.desiredmaps = ["Gen", "Pro", "SC"]
        self.problems = []

        self.getdata()

        self.printwarns(FILES, RunPar)

    def getdata(self):
        file = open(self.full_fname)

        # remove all comments from user (marked with '#', removed by
        # LineFormat). Then, merge back to string to ease map recognition
        lines = []
        for line in file.readlines():
            line = AIM_FH.LineFormat(line)
            if len(line) > 0:
                lines.append(" ".join(line))
        lines = "\n".join(lines)

        # split per map
        temp = lines.split("defmap")
        maps = []
        for mapitem in temp:
            mapitem = mapitem.strip()
            if len(mapitem) > 0:
                maps.append(mapitem)

        for mapitem in maps:
            if mapitem.startswith("Reference"):
                temp = mapitem.split("\n")[1:]
                temp = "\n".join(temp)
                self.references = AIM_RM.readrefstring(temp)
                # temp = temp.split("\n@")
                # for ref in temp:
                #     ref = ref.strip()
                #     if len(ref) > 0:
                #         refparse = AIM_RM.Reference(ref)
                #         if refparse.success:
                #             self.references.append(refparse)
                #     print(refparse)
                #     print(str(refparse))
                #     print("---")
                continue
            if self.maptype == "E":
                submap = E_SubMap(mapitem)
            else:
                submap = D_SubMap(mapitem)
            if submap.success:
                setattr(self, submap.mapname, submap)
            else:
                if submap.mapname in self.desiredmaps:
                    self.success = False
                    self.problems.append([submap.mapname, submap.problem])

    def printwarns(self, FILES, RunPar):
        misstext = "The map is missing."

        # An error should be raised if one of the maps is missing
        for desmap in self.desiredmaps:
            if not hasattr(self, desmap):
                self.success = False
                self.problems.append([desmap, misstext])

        if not self.success:
            ErrorText = (
                "The map stored in the file " + self.full_fname +
                " could not be read due to the following problem(s):\n"
            )
            problems = "\n".join([
                "In defmap " + item[0] + ": " + item[1]
                for item in self.problems
            ])
            fintext = (
                "\nTherefore, this file will be skipped. If you intend to use "
                "this file, this will become an issue soon. In that case, "
                "please make sure this file is formatted correctly."
            )
            ErrorText += problems + fintext

            exhead = " of class AllSourceMaps"

            self.function_warning = AIM_PC.warning(
                ErrorText, False, self.function_warning, FILES.logfilename,
                RunPar, n=3, extra_header=exhead
            )


class E_SubMap:
    def __init__(self, mapitem):

        # Get a list of lines
        rawmap = [item.strip() for item in mapitem.split("\n")]
        tempmap = []
        for item in rawmap:
            if len(item) > 0:
                tempmap.append(item)
        rawmap = tempmap

        self.mapname = rawmap[0].split(" ")[0]

        readmap = read_submap_PEG(rawmap[1:])
        self.success = readmap[0]
        # if the map was not complete/correct
        if not self.success:
            self.problem = readmap[1]
            return
        else:
            self.omega = readmap[1]
            self.PEG = readmap[2]

        self.parse()

    def parse(self):

        # Find out what atoms are relevant
        PEGabs = np.abs(self.PEG)
        relevant_j = []
        for i in range(6):
            if np.amax(PEGabs[i, :]) > 0:
                relevant_j.append(i)
        if len(relevant_j) == 0:
            self.success = False
            self.problem = "The map contains no values, only zeros."
            return
        self.relevant_j = relevant_j
        self.relevant_j_arr = np.array(relevant_j, dtype='int32')

        # Separate out the P, E, and G maps
        self.Pmap = self.PEG[:, 0]
        self.Emap = self.PEG[:, 1:4]
        self.Gmap = self.PEG[:, 4:]

        # See if the gradient is required
        Gabs = np.abs(self.Gmap)
        self.useG = np.amax(Gabs) > 0


class D_SubMap:
    def __init__(self, mapitem):
        self.success = True
        # get a list of lines
        rawmap = [item.strip() for item in mapitem.split("\n")]
        tempmap = []
        for item in rawmap:
            if len(item) > 0:
                tempmap.append(item)
        rawmap = tempmap

        self.mapname = rawmap[0].split(" ")[0]

        len_rawmap = len(rawmap) - 1

        if len_rawmap != 21:
            self.success = False
            self.problem = "The map does not contain the required 21 lines."
            return

        self.omega = []
        self.PEG = []

        for i in range(3):
            readmap = read_submap_PEG(rawmap[(1 + i*7):(1 + (i+1)*7)])
            if not readmap[0]:
                self.success = False
                self.problem = readmap[1]
                return
            else:
                self.omega.append(readmap[1])
                self.PEG.append(readmap[2])

        self.parse()

    def parse(self):

        # create containers to store all directions
        temp_rel_j = set()
        temp_useG = []
        self.Pmap = []
        self.Emap = []
        self.Gmap = []

        for dir in range(3):  # for x, y, z

            # find out what atoms are relevant
            PEGabs = np.abs(self.PEG[dir])

            for i in range(6):
                if np.amax(PEGabs[i, :]) > 0:
                    temp_rel_j.add(i)

            # separate out the P, E, and G maps
            self.Pmap.append(self.PEG[dir][:, 0])
            self.Emap.append(self.PEG[dir][:, 1:4])
            self.Gmap.append(self.PEG[dir][:, 4:])

            Gabs = np.abs(self.Gmap[dir])
            temp_useG.append(np.amax(Gabs) > 0)

        if len(temp_rel_j) == 0:
            self.success = False
            self.problem = "The map contains no values, only zeros."
            return

        self.relevant_j = sorted(list(temp_rel_j))
        self.relevant_j_arr = np.array(self.relevant_j, dtype='int32')

        self.useG = any(temp_useG)


# Here, make classes for each of the coupling maps! These are then called
# by det_Cmap, which also returns them to AllSourceMaps

class NN_Map:
    def __init__(self, full_fname, mapname, FILES, RunPar, function_warning):
        self.full_fname = full_fname
        self.mapname = mapname
        self.function_warning = function_warning
        exhead = " of class AllSourceMaps"
        try:
            self.readmap()
        except Exception:
            self.success = False
        else:
            self.success = True

        # check for completeness:

        desmaps = self.getdesmaps()

        # check if all different maps are present
        if not all(hasattr(self, maptype) for maptype in desmaps):
            self.success = False

        # check if all maps have the right shape/size
        if not all(
            getattr(self, maptype).shape == (13, 13) for maptype in desmaps
        ):
            self.success = False

        if not self.success:
            ErrorText = (
                "The map stored in the file " + self.full_fname +
                " could not be read, as it was of the wrong format, or "
                "incomplete. The file "
                "is expected to contain 21 blocks of 14 lines. The first of "
                "these 14 lines contains only the map name, the 13 after "
                "contain 13 numbers each, separated by spaces. This map comes "
                "with the default AIM installation, and is thus assumed "
                "present and complete. Quitting!"
            )
            self.function_warning = AIM_PC.warning(
                ErrorText, True, self.function_warning, FILES.logfilename,
                RunPar, n=3, extra_header=exhead
            )

    def readmap(self):
        rawmapfile = open(self.full_fname)
        for map_num in range(21):  # there are 21 maps - 3 for each situation!
            mapname = rawmapfile.readline().strip()
            for _ in range(13):  # skip the next 13 lines (dealt with by np)
                _ = rawmapfile.readline()
            mapdata = np.genfromtxt(
                self.full_fname, skip_header=map_num*14 + 1, max_rows=13
            )
            setattr(self, mapname, mapdata)
        rawmapfile.close()

    def getdesmaps(self):
        destypes = ["Coupling", "NtermShift", "CtermShift"]
        desgroups = [""]
        temp = ["Gly_transPro"] + [i+"Pro_transGly" for i in ["", "D"]]
        desgroups.extend([i+j for i in ["_trans", "_cis"] for j in temp])
        desmaps = [i+j for i in destypes for j in desgroups]

        return desmaps


class Tasumi_Map:
    def __init__(self, full_fname, mapname, FILES, RunPar, function_warning):
        self.full_fname = full_fname
        self.mapname = mapname
        self.function_warning = function_warning
        exhead = " of class AllSourceMaps"
        try:
            self.map = np.genfromtxt(self.full_fname, max_rows=13)
        except Exception:
            self.success = False
        else:
            self.success = True

        if not self.map.shape == (13, 13):
            self.success = False

        if not self.success:
            ErrorText = (
                "The map stored in the file " + self.full_fname +
                " could not be read, as it was of the wrong format, or "
                "incomplete. The file is expected to contain 13 lines, which "
                "contain 13 numbers each, separated by spaces. This map comes "
                "with the default AIM installation, and is thus assumed "
                "present and complete. Quitting!"
            )
            self.function_warning = AIM_PC.warning(
                ErrorText, True, self.function_warning, FILES.logfilename,
                RunPar, n=3, extra_header=exhead
            )


class TCC_Map:
    def __init__(self, full_fname, mapname, FILES, RunPar, function_warning):
        self.full_fname = full_fname
        self.mapname = mapname
        self.function_warning = function_warning
        exhead = " of class AllSourceMaps"
        try:
            self.readmap()
        except Exception:
            self.success = False
        else:
            self.success = getattr(self, "success", True)

        if not self.success:
            ErrorText = (
                "The map stored in the file " + self.full_fname +
                " could not be read, as it was of the wrong format, or "
                "incomplete. The file "
                "is expected to contain 18 lines of floats, separated by "
                "spaces. For the exact format, see the file that comes with "
                "AIM upon installation. "
                "This map comes "
                "with the default AIM installation, and is thus assumed "
                "present and complete. Quitting!"
            )
            self.function_warning = AIM_PC.warning(
                ErrorText, True, self.function_warning, FILES.logfilename,
                RunPar, n=3, extra_header=exhead
            )

    def readmap(self):
        rawmap = open(self.full_fname)

        self.fourPiEps = self.looptildata(rawmap)[0]
        self.Gen_alpha, self.Pro_alpha = self.looptildata(rawmap)[:2]

        self.fourPiEps = np.float32(self.fourPiEps)
        self.Gen_alpha = np.float32(self.Gen_alpha)
        self.Pro_alpha = np.float32(self.Pro_alpha)

        self.Gen_q = np.array(self.looptildata(rawmap)[:6], dtype="float32")
        self.Pro_q = np.array(self.looptildata(rawmap)[:6], dtype="float32")
        self.Gen_dq = np.array(self.looptildata(rawmap)[:6], dtype="float32")
        self.Pro_dq = np.array(self.looptildata(rawmap)[:6], dtype="float32")

        temp = []
        for _ in range(6):
            temp.append(self.looptildata(rawmap)[:3])
        self.Gen_v = np.array(temp, dtype="float32")

        temp = []
        for _ in range(6):
            temp.append(self.looptildata(rawmap)[:3])
        self.Pro_v = np.array(temp, dtype="float32")

        for data in ["Gen_q", "Pro_q", "Gen_dq", "Pro_dq"]:
            if getattr(self, data).shape != (6,):
                self.success = False
            else:
                setattr(self, data + "_c", np.ctypeslib.as_ctypes(
                    getattr(self, data)
                ))

        for data in ["Gen_v", "Pro_v"]:
            if getattr(self, data).shape != (6, 3):
                self.success = False
            else:
                setattr(self, data + "_c", AIM_DC.ctype2d(
                    getattr(self, data), 'float32'
                ))

    def looptildata(self, file):
        while True:
            line = file.readline().strip()
            line = line.split('#')[0].strip()
            if len(line) != 0:
                break

        data = [float(i) for i in line.split()]

        return data


class tempcatch:
    def __init__(
        self, full_fname, maptypename, FILES, RunPar, function_warning
    ):
        self.success = False
        self.function_warning = function_warning


class EmptyMap:
    def __init__(self, function_warning):
        self.success = False
        self.function_warning = function_warning


class DmapDud:
    def __init__(self):
        self.Gen = DsubmapDud()
        self.Pro = DsubmapDud()
        self.SC = DsubmapDud()


class DsubmapDud:
    def __init__(self):
        self.useG = False
        self.relevant_j = []


def det_Cmap(full_fname, maptypename, FILES, RunPar, function_warning):

    maptypename = maptypename
    temp = maptypename.split("_")
    # maptype = temp[0]
    mapname = "_".join(temp[1:])

    expected_maps = {
        "NN": NN_Map,
        "Tasumi": Tasumi_Map,
        "TCC": TCC_Map
    }

    if mapname not in expected_maps:
        ErrorText = (
            "Did not recognize the type of the map stored in the file " +
            full_fname +
            ". Only the following names are recognized:\nMap-C_" +
            ".dat\nMap-C_".join(expected_maps.keys()) + ".dat\n" +
            "Therefore, this file will be skipped. If you intended this file "
            "to be one of those listed above, please rename it."
        )
        function_warning = AIM_PC.warning(
            ErrorText, False, function_warning, FILES.logfilename, RunPar,
            n=2, extra_header=" of class AllSourceMaps"
        )
        mapdud = EmptyMap(function_warning)

        return mapdud

    truemap = expected_maps[mapname](
        full_fname, mapname, FILES, RunPar, function_warning
    )

    return truemap


def read_submap_PEG(maplist):
    """
    As both E and D maps contain blocks of 7 lines (1 for omega, 6 for the
    map), here's a single function for both. It takes a list as input. To run
    correctly, this list should be 7 items, one for each line from the file.
    Returns two or three objects.
    1st return: success. True if the submap is correct, False otherwise.
    other returns: if 1st is False, the second is the problem. No third output.
    If 1st is True, 2nd is the value for omega (float), 3rd is PEG (np array)
    """

    # check if the map has the correct length
    if len(maplist) != 7:
        text = "The map does not contain exactly 7 lines."
        return False, text

    # extract omega
    try:
        omega = float(maplist[0])
    except Exception:
        text = "The value for omega is not a float."
        return False, text

    # test if the PEG array is of the correct size
    PEG_list = [item.split() for item in maplist[1:]]
    if any(len(item) != 10 for item in PEG_list):
        text = "Not all lines for the map contain 10 floats"
        return False, text

    # retrieve the PEG grid as list of floats of floats
    try:
        PEG_list = [
            [float(subitem) for subitem in item]
            for item in PEG_list
        ]
    except Exception:
        text = "Not all values for the map are floats"
        return False, text

    # convert the data into a np array
    PEG = np.array(PEG_list, dtype='float32')

    return True, omega, PEG


# ----------------------

class Needed_Maps:
    def __init__(self, RunPar, AllMaps):
        self.Emap = getattr(AllMaps, "E_" + RunPar.map_choice)

        Dmap_choice = RunPar.Dipole_choice
        if Dmap_choice != "Torii":
            self.Dmap = getattr(AllMaps, "D_" + Dmap_choice)
        else:
            self.Dmap = DmapDud()

        for cmap in ["NN", "Tasumi", "TCC"]:
            setattr(self, cmap, getattr(AllMaps, "C_" + cmap))
