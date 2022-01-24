# standard lib imports
import cProfile
import ctypes as ct
import datetime
import os
import sys
import traceback

# 3rd party lib imports
import numpy as np

# my lib imports
import sourcecode.CouplingFunctions as AIM_CF  # just here for .(c)map support
import sourcecode.DataConversion as AIM_DC
import sourcecode.MathFunctions as AIM_MF  # just here for .(c)map support
import sourcecode.PrintCommands as AIM_PC
import sourcecode.PhysicsFunctions as AIM_PF  # just here for .(c)map support
import sourcecode.ReferenceManager as AIM_RM
import sourcecode.SourceMapreader as AIM_SM
import sourcecode.SetOperations as AIM_SO  # just here for .(c)map support
import sourcecode.TimeKeeping as AIM_TK  # just here for .(c)map support
import sourcecode.UniverseMaster as AIM_UM


class FileLocations:
    """
    This class is all about the files AIM reads and writes. The name/location
    of these files are saved as attributes, and its the functions of this class
    that read/write the different files.
    """
    def __init__(self, scriptlocation, script_version, script_version_print):
        """
        Called from AIM.py. Defines the hard coded directory locations and
        crash log file. Cleans out default output and log directories.
        Determines executing OS, and start date/time of the calculation.
        """
        def DirectoryCleanUp(self, nFilesToKeep=30):
            """
            Cleans up the hard-coded log and output directories by removing
            excess files. If the user doesn't change the output location/name,
            these directories will slowly fill up with files otherwise.
            """
            def CleanUp(directory, nFilesToKeep):
                existing_files = [f for f in os.listdir(directory) if
                                  os.path.isfile(os.path.join(directory, f))]
                if len(existing_files) > nFilesToKeep:
                    existing_files.sort()
                    RemoveTheseFiles = existing_files[:-nFilesToKeep]
                    for filename in RemoveTheseFiles:
                        fileloc = os.path.join(directory, filename)
                        os.remove(fileloc)

            CleanUp(self.logdir_hc, nFilesToKeep)
            CleanUp(self.outdir_hc, nFilesToKeep)

        def FindExOS(self):
            """
            Finds out what OS is executing the code, and saves it as
            self.exec_os. The OS information is required when running with a
            not-specified C library: some pre-compiled ones are supplied, but
            compiled c libraries differ per OS...
            """
            # this whole part: defining what OS the script is running on
            if sys.platform == "win32":
                bits = ct.sizeof(ct.c_voidp)*8
                if bits == 32:
                    exec_os = "Windows32bit"
                elif bits == 64:
                    exec_os = "Windows64bit"
                else:
                    ErrorText = (
                        "Environment was determined to be windows, but it is "
                        "neither a 32, nor 64 bit version. It appears to be "
                        + str(bits) + " bit. Please contact the developers to "
                        "solve this. Quitting!")
                    AIM_PC.warning(ErrorText, True, 0, self.logfilename, self)
            elif sys.platform == "darwin":
                exec_os = "MacOS"
            elif sys.platform == "linux":
                exec_os = "Linux"
            else:
                ErrorText = (
                    "executing OS not recognised... sys.platform ="
                    + str(sys.platform)
                    + ". Please contact the developers to solve this. "
                    + "Quitting!")
                AIM_PC.warning(ErrorText, True, 0, self.logfilename, self)
            return exec_os

        # Defining some very basic parameters

        # the full path to the directory in which this script is saved
        self.script_dir = scriptlocation
        self.curr_work_dir = os.getcwd()

        # the current date and time
        self.now = datetime.datetime.now()

        # the current date and time, nicely formatted.
        self.now_str = self.now.strftime("%Y-%m-%d_%H-%M-%S_")

        # the hard-coded log directory
        self.logdir_hc = os.path.join(scriptlocation, "log")

        # the hard-coded output directory
        self.outdir_hc = os.path.join(scriptlocation, "output")

        # the hard-coded sourcefiles directory
        self.sourcedir_hc = os.path.join(scriptlocation, "sourcefiles")

        # the hard-coded extramaps directory
        self.extramapdir_hc = os.path.join(scriptlocation, "extra_maps")

        # a temporary name for the log file, so if the program crashes before
        # it knows where to print the log to, it can print it here, allowing
        # the user to read the error message.
        self.logfilename = os.path.join(
            self.curr_work_dir, self.now_str + "Crash.log")

        # If the program doesn't crash, we don't want to make a crash file and
        # fill it. Therefore, instead of vprinting immediately, we will add any
        # text to this backlog list. Then, if a fatal error occurs, the
        # backlog is printed first, then the error is handled the normal way.
        # if no error occurs, the backlog will be printed to command line and
        # log file as soon as the log file is created.
        # This adding is handled by the vprint command to keep things here
        # universal.
        self.backlog = []

        self.script_version = script_version
        self.script_version_print = script_version_print

        DirectoryCleanUp(self)

        self.exec_os = FindExOS(self)

    def GenerateRunParameters(self, CallCommand):
        """
        Called from AIM.py. Reads the supplied input parameter file, as well
        as the default parameter file (which' location can be changed with the
        input parameter file). From these files, it extracts settings for all
        parameters, and saves these as 'RunPar' - an instance of the
        RunParameters class.
        """

        self.backlog.append([
            1,
            "\n\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\n"
            "/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\\n"
            "\\/\\/\\AmideImaps/\\/\\/\n"
            "/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\\n"
            "\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/"
        ])

        self.backlog.append([
            1,
            "\nRunning version " + self.script_version_print
        ])

        self.backlog.append([
            1,
            "\n" + "*" * 103 +
            "\nThis is an early release of AIM, available while the "
            "corresponding paper is being submitted. The goal\nof this early "
            "release is to aid the scientific community. If you wish to use "
            "this program for your\nscientific work, we'd like to ask you to "
            "be patient with publishing: as soon as the paper is accepted,\n"
            "the first true version will be released. That version will "
            "include a reference and a license.\n" +
            "*" * 103
        ])

        # if a parameter file was supplied, all parameter settings are
        # extracted from it.

        # RawPar will store all parameters as read from the input and defaults
        # files. It will also set the self.Demo_Mode for this FILES class, to
        # use
        RawPar = ParameterFinder(CallCommand, self)

        # The raw parameters are taken and processed: are the choices for
        # settings valid? Do any settings contradict? Do all the requested
        # files actually exist?
        RunPar = RunParameters(RawPar, self)
        if RunPar.use_c_lib:
            self.InitClib()

        # The final used parameters are written to the log file
        self.WriteParToLog(RunPar)
        # An output parameter file is made, to run a future calculation with
        # exactly the same parameters.
        self.WriteParToOutPar(RunPar)

        RunPar.finish_init(self)

        return RunPar

    def set_locations(self, RawPar):
        """
        Called during the RawPar creation, as these attributes should be set
        before RawPar creation completes.
        """
        if hasattr(RawPar, "input_parameter_filename"):
            self.input_parameter_filename = RawPar.input_parameter_filename
        if hasattr(RawPar, "input_parameter_dir"):
            self.input_parameter_dir = RawPar.input_parameter_dir

        # self.Demo_Mode = RawPar.input_parameters.get("Demo_Mode", False)
        self.default_parameter_filename = RawPar.default_parameter_filename
        self.default_parameter_dir = RawPar.default_parameter_dir

        self.sourcedir = RawPar.input_parameters["sourcedir"]
        self.extramapdir = RawPar.input_parameters["extramapdir"]
        self.outdir = RawPar.input_parameters["outdir"]
        self.logdir = RawPar.input_parameters["logdir"]
        self.sourcedir_def = RawPar.sourcedir_def
        self.logfilename = RawPar.input_parameters["logfilename"]
        self.logfilenameraw = RawPar.input_parameters["logfilenameraw"]

        self.Verbose = RawPar.input_parameters.get(
            "Verbose", RawPar.default_parameters["Verbose"]
        )
        self.Verbose_log = RawPar.input_parameters.get(
            "Verbose_log", RawPar.default_parameters["Verbose"]
        )

        for threshold, message in self.backlog:
            AIM_PC.vprint(threshold, message, self.logfilename, self)
        delattr(self, "backlog")

    def FileFinder(self, RawPar, RunPar):
        if RunPar.profiler:
            self.pr = cProfile.Profile()
            self.pr.enable()

        self.outfilenameraw = RawPar.input_parameters["outfilename"]
        self.outfilename = os.path.join(
            self.outdir, NameEdit(self.outfilenameraw, self))
        self.outdipfilenameraw = RawPar.input_parameters["outdipfilename"]
        self.outdipfilename = os.path.join(
            self.outdir, NameEdit(self.outdipfilenameraw, self))
        self.outposfilenameraw = RawPar.input_parameters["outposfilename"]
        self.outposfilename = os.path.join(
            self.outdir, NameEdit(self.outposfilenameraw, self))
        self.outparfilenameraw = RawPar.input_parameters["outparfilename"]
        self.outparfilename = os.path.join(
            self.outdir, NameEdit(self.outparfilenameraw, self))
        self.proffilenameraw = RawPar.input_parameters["proffilename"]
        self.proffilename = os.path.join(
            self.logdir, NameEdit(self.proffilenameraw, self))
        self.pngfilenameraw = RawPar.input_parameters["pngfilename"]
        self.pngfilename = os.path.join(
            self.logdir, NameEdit(self.pngfilenameraw, self))
        self.outramfilenameraw = RawPar.input_parameters["outramfilename"]
        self.outramfilename = os.path.join(
            self.outdir, NameEdit(self.outramfilenameraw, self))

        if "txt" in RunPar.output_format:
            if "Ham" in RunPar.output_type:
                self.outfile = open(self.outfilename + ".txt", "w")
            if "Dip" in RunPar.output_type:
                self.outdipfile = open(self.outdipfilename + ".txt", "w")
            if "Ram" in RunPar.output_type:
                self.outramfile = open(self.outramfilename + ".txt", "w")
            if "Pos" in RunPar.output_type:
                self.outposfile = open(self.outposfilename + ".txt", "w")
        if "bin" in RunPar.output_format:
            if "Ham" in RunPar.output_type:
                self.outfilebin = open(self.outfilename + ".bin", "wb")
            if "Dip" in RunPar.output_type:
                self.outdipfilebin = open(self.outdipfilename + ".bin", "wb")
            if "Ram" in RunPar.output_type:
                self.outramfilebin = open(self.outramfilename + ".bin", "wb")
            if "Pos" in RunPar.output_type:
                self.outposfilebin = open(self.outposfilename + ".bin", "wb")

        # Trying to find functioning map files
        tempdict = {}
        for item in [
            "referencefile", "resnamesfile", "atnamesfile"
        ]:
            if item not in RawPar.input_parameters:
                tempdict[item+"filename"] = os.path.join(
                    self.sourcedir_def, RawPar.default_parameters[item])
                try:
                    _ = open(tempdict[item+"filename"])
                except FileNotFoundError:
                    tempdict[item+"filename"] = os.path.join(
                        self.sourcedir_hc, RawPar.default_parameters[item])
                    try:
                        _ = open(tempdict[item+"filename"])
                    except FileNotFoundError:
                        ErrorText = (
                            "\nThe following sourcefile could not be found: "
                            + item + ".\nPlease refer to the manual, as this "
                            + "file is required for calculations. Quitting!")
                        AIM_PC.warning(
                            ErrorText, True, 0, self.logfilename, RunPar)
            else:
                tempdict[item+"filename"] = os.path.join(
                    self.sourcedir, RawPar.input_parameters[item])
                try:
                    _ = open(tempdict[item+"filename"])
                except FileNotFoundError:
                    ErrorText = (
                        "\nThe following sourcefile could not be found: "
                        + item + ".\nThe input file specifies an alternative "
                        + "location for this file. Either make sure it is "
                        + "present in this new location, or disable the line "
                        + "from the input file to use the default file. "
                        + "Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, 0, self.logfilename, RunPar)

        # self.Tokmakoff_EMapfilename = tempdict["Tokmakoff_EMapfilename"]
        # self.Skinner_EMapfilename = tempdict["Skinner_EMapfilename"]
        # self.Jansen_EMapfilename = tempdict["Jansen_EMapfilename"]
        # self.Custom_EMapfilename = tempdict["Custom_EMapfilename"]
        # self.Jansen_DMapfilename = tempdict["Jansen_DMapfilename"]
        # self.NN_Mapfilename = tempdict["NN_Mapfilename"]
        # self.Tasumi_Mapfilename = tempdict["Tasumi_Mapfilename"]
        # self.TCC_Mapfilename = tempdict["TCC_Mapfilename"]
        self.referencefilefilename = tempdict["referencefilefilename"]
        self.resnamesfilefilename = tempdict["resnamesfilefilename"]
        self.atnamesfilefilename = tempdict["atnamesfilefilename"]

        # for finding top and trj files. If neither are specified, enter demo
        # mode. If one is specified, thats a problem. If both are specified,
        # use them!
        temp = 0
        for item in ["topfile", "trjfile"]:
            if item in RawPar.input_parameters:
                temp += 1

        if temp == 1:
            for item in ["topfile", "trjfile"]:
                if item not in RawPar.input_parameters:
                    ErrorText = (
                        "\nInput file is missing the following required "
                        + "parameter: " + item + "\n Either specify it, or "
                        + "make sure to omit both 'topfile' and 'trjfile' "
                        + "parameters to run in Demo Mode. Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, 0, self.logfilename, RunPar)
        elif temp == 2:
            self.topfile = os.path.join(
                self.input_parameter_dir, RawPar.input_parameters["topfile"])
            self.trjfile = os.path.join(
                self.input_parameter_dir, RawPar.input_parameters["trjfile"])
        elif temp == 0:
            self.topfile = os.path.join(
                self.default_parameter_dir,
                RawPar.default_parameters["topfile"])
            self.trjfile = os.path.join(
                self.default_parameter_dir,
                RawPar.default_parameters["trjfile"])
            self.Demo_Mode = True

    def ClibFinder(self, RawPar, RunPar):
        def TryLoadUnspecCLib(self, RawPar):
            returnval = True
            libfilename = os.path.join(
                self.sourcedir_def, RawPar.default_parameters["libfile"])
            try:
                _ = ct.CDLL(libfilename)
            except Exception:
                if self.exec_os == "Windows32bit":
                    libfileOS = RawPar.default_parameters["libfile_Win32"]
                elif self.exec_os == "Windows64bit":
                    libfileOS = RawPar.default_parameters["libfile_Win64"]
                elif self.exec_os == "MacOS":
                    libfileOS = RawPar.default_parameters["libfile_MacOS"]
                elif self.exec_os == "Linux":
                    libfileOS = RawPar.default_parameters["libfile_Linux"]
                libfilename = os.path.join(self.sourcedir_def, libfileOS)
                try:
                    _ = ct.CDLL(libfilename)
                except Exception:
                    try:
                        libfilename = os.path.join(
                            self.sourcedir_hc, libfileOS)
                        _ = ct.CDLL(libfilename)
                    except Exception:
                        returnval = False
            return returnval, libfilename

        if "use_c_lib" not in RawPar.input_parameters:
            if "libfile" not in RawPar.input_parameters:
                success, libfilename = TryLoadUnspecCLib(self, RawPar)
                if success:
                    ErrorText = (
                        "\nInput file did not specify whether to use a "
                        "c-library for speedup, but a functioning one was "
                        "found and will thus be used!")
                    AIM_PC.vprint(1, ErrorText, self.logfilename, RunPar)
                else:
                    ErrorText = (
                        "\nNo working c-library file found. If you would like "
                        "to use the c-lib speedup, see the manual. Running "
                        "without c-lib speedup!")
                    AIM_PC.vprint(1, ErrorText, self.logfilename, RunPar)
            else:
                libfilename = os.path.join(
                    self.sourcedir, RawPar.input_parameters["libfile"])
                try:
                    _ = ct.CDLL(libfilename)
                except Exception:
                    ErrorText = (
                        "\nThe specified c-library file was not found. If you "
                        "would like to use the c-lib speedup, see the manual. "
                        "Running without c-lib speedup!")
                    AIM_PC.vprint(1, ErrorText, self.logfilename, RunPar)
                    success = False
                else:
                    ErrorText = (
                        "\nThe specified c-library file has been found and is "
                        "functional. As the input parameter file did not "
                        "specify whether to use it, and it works, it is used "
                        "anyways.")
                    AIM_PC.vprint(1, ErrorText, self.logfilename, RunPar)
                    success = True
        elif RawPar.input_parameters["use_c_lib"]:
            if "libfile" not in RawPar.input_parameters:
                success, libfilename = TryLoadUnspecCLib(self, RawPar)
                if not success:
                    ErrorText = (
                        "\nNo working c-library file found. As the use of one "
                        "was specifically requested, make sure one is "
                        "available for your system. Instructions on how to do "
                        "so can be found in the manual. Alternatively, set "
                        "'use_c_lib' to False, or don't specify it at all for "
                        "a non-c-speedup run. Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, 0, self.logfilename, RunPar)
            else:
                libfilename = os.path.join(
                    self.sourcedir, RawPar.input_parameters["libfile"])
                try:
                    _ = ct.CDLL(libfilename)
                except Exception:
                    ErrorText = (
                        "\nThe specified c-library file was not found. As the "
                        "use of one was specifically requested, make sure one "
                        "is available for your system. Instructions on how to "
                        "do so can be found in the manual. Alternatively, "
                        "don't specify 'libfile' in the input parameter file, "
                        "set 'use_c_lib' to False, or don't specify it at all "
                        "for a non-c-speedup run. Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, 0, self.logfilename, RunPar)
                else:
                    success = True
        else:
            success = False

        if success:
            self.libfile = libfilename
            self.clib = ct.CDLL(libfilename)

        return success

    def InitClib(self):

        # CalcFieldGrad
        self.clib.CalcFieldGrad.argtypes = [
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.c_int, ct.POINTER(ct.c_float)
        ]
        self.clib.CalcFieldGrad.restype = None

        # CalcField
        self.clib.CalcField.argtypes = [
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.c_int, ct.POINTER(ct.c_float)
        ]
        self.clib.CalcField.restype = None

        # group_difference
        self.clib.group_difference.argtypes = [
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)
        ]
        self.clib.group_difference.restype = None

        # res_finder
        self.clib.res_finder.argtypes = [
            ct.POINTER(ct.c_int), ct.c_int, ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.res_finder.restype = None

        # res_finder_dumb
        self.clib.res_finder_dumb.argtypes = [
            ct.POINTER(ct.c_int), ct.c_int, ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.res_finder_dumb.restype = None

        # AGsorter
        self.clib.AGsorter.argtypes = [
            ct.c_int, ct.POINTER(ct.c_int), ct.c_int, ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_int), ct.c_float, ct.c_int, ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float)
        ]
        self.clib.AGsorter.restype = np.int32

        # ix_builder
        self.clib.ix_builder.argtypes = [
            ct.POINTER(ct.c_int), ct.c_int, ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)
        ]
        self.clib.ix_builder.restype = np.int32

        # PrepTCC
        self.clib.PrepTCC.argtypes = [
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_int), ct.c_int,
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.c_float, ct.c_float, ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.PrepTCC.restype = None

        # CalcTCC
        self.clib.CalcTCC.argtypes = [
            ct.c_int, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.c_float, ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.CalcTCC.restype = ct.c_float

        # PrepTDCKrimm
        self.clib.PrepTDCKrimm.argtypes = [
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_int), ct.c_int,
            ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.PrepTDCKrimm.restype = None

        # CalcTDCKrimm
        self.clib.CalcTDCKrimm.argtypes = [
            ct.c_int, ct.c_int, ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.CalcTDCKrimm.restype = ct.c_float

        # PrepTDCTasumi
        self.clib.PrepTDCTasumi.argtypes = [
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_int), ct.c_int,
            ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.PrepTDCTasumi.restype = None

        # CalcTDCTasumi
        self.clib.CalcTDCTasumi.argtypes = [
            ct.c_int, ct.c_int, ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.CalcTDCTasumi.restype = ct.c_float

        # CalcGenCoup
        self.clib.CalcGenCoup.argtypes = [
            ct.c_int, ct.c_int, ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
        ]
        self.clib.CalcGenCoup.restype = ct.c_float

    def WriteParToLog(self, RunPar):
        logfile = open(self.logfilename, "a")

        logfile.write(
            "\n\n\n====================\nFinal input "
            "parameters\n====================\n")
        logfile.write("\nUsed files:")

        if hasattr(self, "input_parameter_filename"):
            logwrite(
                logfile, "Input parameter file", 3,
                self.input_parameter_filename)
        logwrite(logfile, "Input topology file", 3, self.topfile)
        logwrite(logfile, "Input trajectory file", 3, self.trjfile)
        if "Ham" in RunPar.output_type:
            logwrite(logfile, "Output hamiltonian file", 3, self.outfilename)
        else:
            logfile.write("\nNo Hamiltonian calculated")

        if "Dip" in RunPar.output_type:
            logwrite(logfile, "Output dipole file", 3, self.outdipfilename)
        else:
            logfile.write("\nNo dipoles calculated")

        logwrite(logfile, "Output parameter file", 3, self.outparfilename)

        if "Pos" in RunPar.output_type:
            logwrite(logfile, "Output position file", 3, self.outposfilename)
        else:
            logfile.write("\nNo positions saved")

        if "Ram" in RunPar.output_type:
            logwrite(logfile, "Output raman file", 3, self.outramfilename)
        else:
            logfile.write("\nNo Raman tensor calculated")

        if RunPar.use_c_lib:
            logfile.write(
                "\n\n   For further speedup, this run utilizes the following "
                + "c++ library:\n   " + self.libfile)
        else:
            logfile.write("\n\n   Using numba speedup")

        if RunPar.profiler:
            if RunPar.pngout:
                logfile.write(
                    "\n\n   This run is profiled for further runtime "
                    "information, and a graphical representation will be "
                    "saved at:")
            else:
                logfile.write(
                    "\n\n   This run is profiled for further runtime "
                    "information, but no graphical representation will be "
                    "saved")
            logwrite(logfile, "   Profiler raw output file", 3,
                     self.proffilename)
            if RunPar.pngout:
                logwrite(logfile, "   Profiler png output file", 3,
                         self.pngfilename)
        else:
            logfile.write(
                "\n   This run is not being profiled. Start and end times will"
                " still be recorded")

        logfile.write("\n\n\nRun parameters:")
        logwrite(logfile, "Influencers selection", 3, "")
        for item in RunPar.influencers:
            logfile.write(item + " ")
        logwrite(logfile, "Oscillators selection", 3, "")
        for item in RunPar.oscillators:
            logfile.write(item + " ")
        logwrite(logfile, "Applying dipole-dipole coupling to",
                 3, RunPar.apply_dd_coupling)
        logwrite(logfile, "Electrostatic map", 3, RunPar.map_choice)
        logwrite(logfile, "Dipole method", 3, RunPar.Dipole_choice)
        logwrite(logfile, "Coupling method", 3, RunPar.coupling_choice)
        logwrite(logfile, "NN coupling method", 3, RunPar.NN_coupling_choice)
        logwrite(logfile, "Positions of atom", 3, RunPar.AtomPos_choice)
        logwrite(logfile, "Starting frame number", 3, str(RunPar.start_frame))
        logwrite(logfile, "Amount of frames to calculate", 3,
                 str(RunPar.nFrames_to_calculate))
        logwrite(logfile, "End frame number", 3, str(RunPar.end_frame))
        logwrite(logfile, "Cutoff distance for Electrostatic effects", 3,
                 str(RunPar.SphereSize) + " Angstrom\n")

        if RunPar.replicate_orig_AIM:
            logfile.write("\n   Replicating behaviour of "
                          "original AmideImaps.c script")
        else:
            logfile.write("\n   Not replicating behaviour of original "
                          "AmideImaps.c script")

        if RunPar.NSA_toggle:
            logfile.write(
                "\n   Using a faster method of finding atoms within cut-off "
                "distance, with the following parameters:")
            logwrite(
                logfile, "Cutoff distance for buffer group", 6,
                str(RunPar.NSA_spheresize) + " Angstrom")
            logwrite(
                logfile, "Buffer group lifetime", 6, str(RunPar.NSA_nframes)
                + " Frames \n")
        else:
            logfile.write(
                "\n   No faster method for finding atoms within cut-off "
                "distance used")

        if RunPar.atom_based_chainID:
            logfile.write(
                "\n   Identifying protein chains using atoms and bonds"
            )
        else:
            logfile.write(
                "\n   Identifying protein chains using residue, molecule and "
                "segment numbers in topology file"
            )

        if RunPar.use_protein_specials:
            logfile.write(
                "\n   Considering the special aminoacid names a part of the "
                "protein")
        else:
            logfile.write(
                "\n   Only allowing the classic 20 aminoacids to make up the "
                "protein")

        logfile.write("\n   All calculated long-range coupling constants "
                      "are multiplied by " + str(RunPar.Scale_LRCoupling))

        if RunPar.Use_AmGroup_selection_criteria:
            logfile.write(
                "\n\n   Applying further restrictions on which amide groups "
                "to include in the hamiltonian and dipoles")
            if RunPar.Use_resname_blacklist:
                logwrite(logfile, "Residue name blacklist", 6, "")
                for item in RunPar.resname_blacklist:
                    logfile.write(item + " ")
            if RunPar.Use_resname_whitelist:
                logwrite(logfile, "Residue name whitelist", 6, "")
                for item in RunPar.resname_whitelist:
                    logfile.write(item + " ")
            if RunPar.Use_resnum_blacklist:
                logwrite(logfile, "Residue number blacklist", 6, "")
                for item in RunPar.resnum_blacklist:
                    logfile.write(str(item) + " ")
            if RunPar.Use_resnum_whitelist:
                logwrite(logfile, "Residue number whitelist", 6, "")
                for item in RunPar.resnum_whitelist:
                    logfile.write(str(item) + " ")
        else:
            logfile.write(
                "\n   Not applying any further restrictions on which amide "
                "groups to include in the hamiltonian and dipoles")

        if RunPar.TreatNN:
            logfile.write("\n   Considering Nearest Neigbour effects")
        else:
            logfile.write("\n   Not considering Nearest Neigbour effects")

        logfile.write("\n\n")

        logfile.close()

    def WriteParToOutPar(self, RunPar):
        outparfile = open(self.outparfilename, "w")

        outparfile.write("# topology and trajectory files, including folder\n")
        bulkparwrite(outparfile, 0, ["topfile", "trjfile"], self)

        outparfile.write("\n\n\n# Support files location\n")
        templist = [
            "sourcedir", "referencefilefilename", "resnamesfilefilename",
            "atnamesfilefilename"]
        if RunPar.use_c_lib:
            templist.append("libfile")

        for index, value in enumerate(templist):
            if index == 0:
                parwrite(outparfile, value, 0, getattr(self, value))
                continue
            if len(value) < 8 or value[-8:] != "filename":
                temp = value
            else:
                temp = value[:-8]
            try:
                temp2 = os.path.relpath(getattr(self, value),
                                        start=self.sourcedir)
            except Exception:
                parwrite(outparfile, temp, 0, getattr(self, value))
            else:
                parwrite(outparfile, temp, 0, temp2)

        parwrite(outparfile, "use_c_lib", 0, RunPar.use_c_lib)

        parwrite(outparfile, "extramapdir", 0, self.extramapdir)

        outparfile.write("\n\n\n\n# output files\n")
        bulkparwrite(outparfile, 0, [
            "outdir", "outfilenameraw", "outdipfilenameraw",
            "outposfilenameraw", "outramfilenameraw",  "outparfilenameraw"
        ], self)
        outparfile.write("\n")
        bulkparwrite(outparfile, 0, [
            "logdir", "logfilenameraw", "proffilenameraw",
            "pngfilenameraw"
        ], self)

        outparfile.write("\n\n\n\n# Output file settings\n")
        bulkparwrite(outparfile, 0, [
            "output_format", "output_type", "Verbose", "Verbose_log",
            "profiler", "pngout"
        ], RunPar)

        outparfile.write("\n\n\n\n# Run Parameters\n")
        bulkparwrite(outparfile, 0, [
            "influencers", "oscillators", "apply_dd_coupling", "map_choice",
            "Dipole_choice", "coupling_choice",
            "NN_coupling_choice", "AtomPos_choice", "start_frame",
            "nFrames_to_calculate", "end_frame", "max_time", "SphereSize"
        ], RunPar)
        outparfile.write("\n")
        bulkparwrite(outparfile, 0, [
            "replicate_orig_AIM", "NSA_toggle", "NSA_nframes",
            "NSA_spheresize"
        ], RunPar)
        outparfile.write("\n")
        bulkparwrite(outparfile, 0, [
            "atom_based_chainID", "use_protein_specials", "Scale_LRCoupling"
        ], RunPar)
        outparfile.write("\n")
        bulkparwrite(outparfile, 0, [
            "Use_AmGroup_selection_criteria", "resnum_whitelist",
            "resnum_blacklist", "resname_whitelist", "resname_blacklist"
        ], RunPar)
        outparfile.write("\n")
        bulkparwrite(outparfile, 0, ["TreatNN"], RunPar)

        outparfile.close()

    def InitProgram(self, RunPar):
        if RunPar.profiler:
            AIM_PC.vprint(
                2, "\n\nThis run is being profiled! This slows down the "
                "program considerably", self.logfilename, RunPar)

        # Reads in the entire system. The resulting WS is independent of the MD
        # package used: all conflicts are solved here. All basic properties are
        # read (like charge, mass, position) and converted to a c-friendly
        # format if necessary.
        WS = AIM_UM.Universe(self, RunPar)

        return WS

    def ReadMaps(self, RunPar):
        self.AllMaps = AIM_SM.AllSourceMaps(self, RunPar)
        # AllMaps = MapReader_old(self, RunPar)
        # return AllMaps
        self.OtherRefs = AIM_RM.readreffile(self)

    def WriteOutput(self, RunPar, WS):
        if "bin" in RunPar.output_format:  # binary output
            framenum_arr = np.array([WS.framenum], dtype='float32')

            if "Ham" in RunPar.output_type:
                Hamil_bin = np.empty(
                    (int(WS.res_desired_len*(WS.res_desired_len+1)/2)),
                    dtype='float32')
                temp = 0
                for i in range(WS.res_desired_len):
                    for j in range(i, WS.res_desired_len):
                        Hamil_bin[temp] = WS.Hamiltonian[i, j]
                        temp += 1
                framenum_arr.tofile(self.outfilebin)
                Hamil_bin.tofile(self.outfilebin)

            if "Dip" in RunPar.output_type:
                Dip_bin = np.empty((WS.res_desired_len*3), dtype='float32')
                for i in range(3):
                    for j in range(WS.res_desired_len):
                        Dip_bin[i*WS.res_desired_len + j] = WS.AtomPos[j, i]
                framenum_arr.tofile(self.outdipfilebin)
                Dip_bin.tofile(self.outdipfilebin)

            if "Ram" in RunPar.output_type:
                Ram_bin = np.empty((WS.res_desired_len*6), dtype='float32')
                for i in range(6):
                    for j in range(WS.res_desired_len):
                        Ram_bin[i*WS.res_desired_len + j] = WS.Raman[j, i]
                framenum_arr.tofile(self.outramfilebin)
                Ram_bin.tofile(self.outramfilebin)

            if "Pos" in RunPar.output_type:
                Pos_bin = np.empty((WS.res_desired_len*3), dtype='float32')
                for i in range(3):
                    for j in range(WS.res_desired_len):
                        Pos_bin[i*WS.res_desired_len + j] = WS.AtomPos[j, i]
                framenum_arr.tofile(self.outposfilebin)
                Pos_bin.tofile(self.outposfilebin)

        if "txt" in RunPar.output_format:  # text output
            # rounding all values for writing them to file.
            if "Ham" in RunPar.output_type:
                Hamiltonian = np.round(WS.Hamiltonian, decimals=6)
                self.outfile.write(str(WS.framenum) + " ")
                for i in range(WS.res_desired_len):
                    for j in range(i, WS.res_desired_len):
                        self.outfile.write(str(Hamiltonian[i, j]) + " ")
                self.outfile.write("\n")

            if "Dip" in RunPar.output_type:
                Dipoles = np.round(WS.Dipoles, decimals=6)
                self.outdipfile.write(str(WS.framenum) + " ")
                for i in range(3):
                    for j in range(WS.res_desired_len):
                        self.outdipfile.write(str(Dipoles[j, i]) + " ")
                self.outdipfile.write("\n")

            if "Ram" in RunPar.output_type:
                Raman = np.round(WS.Raman, decimals=6)
                self.outramfile.write(str(WS.framenum) + " ")
                for i in range(6):
                    for j in range(WS.res_desired_len):
                        self.outramfile.write(str(Raman[j, i]) + " ")
                self.outramfile.write("\n")

            if "Pos" in RunPar.output_type:
                AtomPos = np.round(WS.AtomPos, decimals=6)
                self.outposfile.write(str(WS.framenum) + " ")
                for i in range(3):
                    for j in range(WS.res_desired_len):
                        self.outposfile.write(str(AtomPos[j, i]) + " ")
                self.outposfile.write("\n")


class ParameterFinder:
    """
    This class deals with the initial choices of parameters. It reads them from
    the files and checks for completeness. The final choices are saved in a
    RunParameters class instance.
    """
    def __init__(self, CallCommand, FILES):
        """
        Called from FileLocations.GenerateRunParameters. Calls all functions
        of this class, as well as a single one from the FileLocations class.
        """

        self.more_inp_setnames = [
            "outdir", "logdir"
        ]
        self.ReadInputs(CallCommand, FILES)

        # read the default parameter file. It is possible another default
        # was specified in the input file, this is taken into account.
        # All settings are read from the default file, and some are immediately
        # taken from RawPar into FILES.
        setattr(
            FILES, "Demo_Mode",
            self.input_parameters.get("Demo_Mode", False)
        )
        self.more_def_setnames = [
            "libfile_Win32", "libfile_Win64",
            "libfile_MacOS", "libfile_Linux"
        ]
        self.ReadDefaults(FILES)

        # Check if the default parameter file is complete (all var's present)
        self.CheckDefaults(FILES)

        # See if the requested log/output/sourcefiles folders exist.
        self.Locationcheck(FILES)

        # Setting some attributes of FILES before we can continue
        FILES.set_locations(self)

        # The specified parameters in the input file are written to the log
        # file, before solving any issues.
        self.WriteRawInputToLog(FILES)
        # Any parameters missing in the input file are supplied from the
        # default file.
        self.ParameterCheck(FILES)

    def ReadInputs(self, CallCommand, FILES):
        """
        Based on the user input in the commandline, it determines the input
        parameter file. If an input parameter file was supplied, it reads the
        contents and saves those as an attribute of itself.
        """
        if len(CallCommand) > 1:
            # define some useful parameters:

            # name of input parameter file, including full path
            self.input_parameter_filename = os.path.realpath(CallCommand[1])

            # directory in which the input parameter file is saved
            self.input_parameter_dir = os.path.dirname(
                self.input_parameter_filename)

            # extract the parameters and their settings from the parameter file
            (self.input_parameters,
             self.input_parameters_str,
             self.unrecognised_parameters) = read_parameter_file(
                 self.input_parameter_filename, FILES, "input",
                 more_setname_str=self.more_inp_setnames)
        else:
            FILES.backlog.append([
                0,
                "\nNo input parameter file given. Entering Demo Mode!\n"
            ])
            self.input_parameters = {"Demo_Mode": True}
            self.input_parameters_str, self.unrecognised_parameters = {}, []

    def ReadDefaults(self, FILES):
        # for finding the default parameter file. If a name was specified in
        # input file, use that. Use hard-coded default otherwise.
        if "def_parfile" not in self.input_parameters:
            # no need to check whether sourcedir exists, its already within a
            # try anyways.

            try:
                temp = os.path.join(
                    self.input_parameter_dir,
                    self.input_parameters["sourcedir"],
                    "default_input.txt")
                _ = open(temp)
            except Exception:
                temp = os.path.join(FILES.sourcedir_hc, "default_input.txt")
                try:
                    _ = open(temp)
                except FileNotFoundError:
                    ErrorText = ("The default parameter file cannot be found. "
                                 "Please refer to the manual, as this file is "
                                 "required for calculations. Quitting!")
                    AIM_PC.warning
                    (ErrorText, True, 0, FILES.logfilename, FILES)
                else:
                    self.input_parameters["def_parfilename"] = temp
            else:
                self.input_parameters["def_parfilename"] = temp
        else:
            if "sourcedir" in self.input_parameters:
                temp = os.path.join(
                    self.input_parameter_dir,
                    self.input_parameters["sourcedir"],
                    self.input_parameters["def_parfile"])
            else:
                temp = os.path.join(
                    FILES.sourcedir_hc,
                    self.input_parameters["def_parfile"])
            try:
                _ = open(temp)
            except FileNotFoundError:
                ErrorText = (
                    "The default parameter file cannot be found. The input "
                    "file specifies an alternative location for this file. "
                    "Either make sure it is present in this new location, or "
                    "disable the line from the input file to use the default "
                    "file. Quitting!")
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, FILES)
            else:
                self.input_parameters["def_parfilename"] = temp

        # now the file has been found, define location.

        # the whole path + name
        self.default_parameter_filename = os.path.realpath(
            self.input_parameters["def_parfilename"])
        # the full path to the directory in which the file lives
        self.default_parameter_dir = os.path.dirname(
            self.default_parameter_filename)

        # read-in default parameters
        self.default_parameters = read_parameter_file(
            self.default_parameter_filename, FILES, "default",
            more_setname_str=self.more_def_setnames
            )[0]

    def CheckDefaults(self, FILES):
        # check if the default parameter file is complete
        SetNames = get_settings_names()
        all_settings = (
            SetNames["settings_names_list_all"] + self.more_def_setnames
        )
        for item in all_settings:
            if item not in self.default_parameters:
                if (item in ["NSA_spheresize", "NSA_nframes"]
                        and not self.default_parameters["NSA_toggle"]):
                    continue
                if (item in SetNames["settings_names_list_BWlist"]
                        or item in ["def_parfile"]):
                    continue
                ErrorText = (
                    "The following parameter is missing from default "
                    "parameters file: " + item + "\nThe default parameters "
                    "file should contain all parameters. Please provide all, "
                    "and retry. Quitting!")
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, FILES)

    def Locationcheck(self, FILES):
        # These are no longer expected in default_parameters, as the user
        # must give them in the input_parameters file. If not, they should
        # default to the folder from which the program is run.
        for item in ["outdir", "logdir"]:
            self.default_parameters[item] = FILES.curr_work_dir

        for item in ["sourcedir", "extramapdir", "outdir", "logdir"]:
            if item not in self.input_parameters:
                self.input_parameters[item] = os.path.join(
                    self.default_parameter_dir, self.default_parameters[item])
            else:
                self.input_parameters[item] = os.path.join(
                    self.input_parameter_dir, self.input_parameters[item])
            self.input_parameters[item] = os.path.realpath(
                self.input_parameters[item])
        self.sourcedir_def = os.path.join(
            self.default_parameter_dir, self.default_parameters["sourcedir"])

        if "logfilename" in self.input_parameters:
            # fetch from input
            self.input_parameters["logfilenameraw"] = (
                self.input_parameters["logfilename"])
        else:
            # fetch from default
            self.input_parameters["logfilenameraw"] = (
                self.default_parameters["logfilename"])
        temp = NameEdit(self.input_parameters["logfilenameraw"], FILES)
        self.input_parameters["logfilename"] = os.path.join(
            self.input_parameters["logdir"], temp)

    def WriteRawInputToLog(self, FILES):
        # start writing the log file
        try:
            logfile = open(FILES.logfilename, "a")
        except Exception:
            logfile = open(FILES.logfilename, "w")
        logfile.write(
            "====================\nInput file\n====================\n")

        # write the raw inputfile to the log file
        for key, value in self.input_parameters_str.items():
            logwrite(logfile, key, 0, value)

        logfile.close()

    def ParameterCheck(self, FILES):
        logfile = open(FILES.logfilename, "a")
        function_warning = 0
        logfile.write(
            "\n\n\n====================\nChecking "
            "parameters\n====================\n")

        # now, make sure all input parameters make sense/are complete (do a
        # proper warning if not, so this can be seen in log file). Then, print
        # the final set of parameters to log file in the same way as previous
        # versions, and write an input-safe version to a new file in the
        # outfolder. This new file MUST be directly usable as input file for
        # new calculation.

        # if the input file contained any unknown parameter names, an extra
        # warning is printed to the log file. Although a message was printed
        # to the command line (and temporary log file?) for each of these
        # parameters, this can be missed by the user. This makes sure that
        # this info can be retrieved later from log file as well!

        if len(self.unrecognised_parameters) != 0:
            logfile.write(
                "\n\nUh-oh! Not all parameter names from the input file were "
                "recognised! The following were unclear to the program:")
            for parname in self.unrecognised_parameters:
                logfile.write(
                    "\n" + parname + self.input_parameters_str[parname])

        # The code can run without a specified setting for these parameters
        # without any issues. #Raman#
        for item in [
            "outfilename", "outdipfilename", "outramfilename",
            "outposfilename",
            "outparfilename", "proffilename", "pngfilename",
            "output_format", "output_type", "Verbose",
            "Verbose_log", "profiler", "pngout", "influencers", "oscillators",
            "apply_dd_coupling", "AtomPos_choice", "max_time",
            "replicate_orig_AIM", "NSA_toggle", "atom_based_chainID",
            "use_protein_specials", "Scale_LRCoupling", "TreatNN"
        ]:
            if item not in self.input_parameters:
                self.input_parameters[item] = self.default_parameters[item]

        logfile.close()
        # Although the code is technically able to run without these parameters
        # specified, they are too important not to notify the user they're
        # missing.
        for item in [
            "map_choice", "Dipole_choice", "coupling_choice",
            "NN_coupling_choice", "SphereSize"
        ]:
            if item not in self.input_parameters:
                self.input_parameters[item] = self.default_parameters[item]
                ErrorText = (
                    "\nInput file is missing the following parameter: "
                    + item + "\nThe following default setting will be used: "
                    + str(self.default_parameters[item]) +
                    "\nBe warned that this might not be suitable to your "
                    + "system!\n")
                if not FILES.Demo_Mode:
                    function_warning = AIM_PC.warning(
                        ErrorText, False, function_warning, FILES.logfilename,
                        FILES)

        AIM_PC.finwarning(function_warning, FILES.logfilename, FILES)


class RunParameters:
    """
    This class contains all the final choices of parameters. It takes the input
    choices from a ParameterFinder class instance, and parses the choices. It
    checks whether choices for different parameters are valid, and fixes any
    possible contradictions between parameters.
    """
    def __init__(self, RawPar, FILES):
        # migrate settings into RunPar
        self.InitRunSet(RawPar, FILES)
        # check if all specified files exist
        FILES.FileFinder(RawPar, self)
        self.Demo_Mode = FILES.Demo_Mode

        # check if the specified c-library exists
        success = FILES.ClibFinder(RawPar, self)
        self.use_c_lib = success
        if self.use_c_lib:
            self.SphereSize_c = np.float32(self.SphereSize)

        # interpret the black-/whitelist settings
        self.InitBWlist(RawPar, FILES.logfilename)
        # interpret the start/end frame settings
        self.InitFrames(RawPar, FILES.logfilename)
        # make sure the settings for the NSA make sense
        self.InitNSA(RawPar, FILES.logfilename)
        # fix contradicting settings (like requesting to not calculate
        # anything)
        self.FixContradictions(FILES.logfilename, FILES)
        # see if the choices for settings are valid (are map-name choices
        # recognized? Are booleans selected using true/false? etc.)
        self.ValidateInput(FILES.logfilename)

        # read in the data for the source-maps, and validate
        FILES.ReadMaps(self)
        self.ValidateMaps(FILES)

        # read in the data for (c)maps
        # However, after all available maps have been read, the program must
        # still check whether the choice of oscillating groups is valid in
        # the input parameter file! Don't forget to do so here!

        self.ReadExtraMaps(FILES)
        self.ValidateExtraMaps(FILES)

    def InitRunSet(self, RawPar, FILES):

        self.output_format = RawPar.input_parameters["output_format"]
        self.output_type = RawPar.input_parameters["output_type"]

        self.Verbose = RawPar.input_parameters["Verbose"]
        self.Verbose_log = RawPar.input_parameters["Verbose_log"]

        self.profiler = RawPar.input_parameters["profiler"]
        if self.profiler:
            self.pngout = RawPar.input_parameters["pngout"]
        else:
            self.pngout = False
            if RawPar.input_parameters["pngout"]:
                logfile = open(FILES.logfilename, "a")
                logfile.write(
                    "\n\nAs the run will not be profiled, 'pngout' is set to "
                    + "'false': profiling is necessary for generating a "
                    + "runtime graph.")

        self.influencers = RawPar.input_parameters["influencers"]
        self.oscillators = RawPar.input_parameters["oscillators"]
        self.apply_dd_coupling = RawPar.input_parameters["apply_dd_coupling"]
        self.map_choice = RawPar.input_parameters["map_choice"]
        self.Dipole_choice = RawPar.input_parameters["Dipole_choice"]
        self.coupling_choice = RawPar.input_parameters["coupling_choice"]
        self.NN_coupling_choice = RawPar.input_parameters["NN_coupling_choice"]
        self.AtomPos_choice = RawPar.input_parameters["AtomPos_choice"]

        self.max_time = RawPar.input_parameters["max_time"]
        self.SphereSize = RawPar.input_parameters["SphereSize"]

        self.replicate_orig_AIM = RawPar.input_parameters["replicate_orig_AIM"]

        self.atom_based_chainID = RawPar.input_parameters["atom_based_chainID"]
        self.use_protein_specials = (
            RawPar.input_parameters["use_protein_specials"])
        self.Scale_LRCoupling = RawPar.input_parameters["Scale_LRCoupling"]

        self.TreatNN = RawPar.input_parameters["TreatNN"]

        self.Demo_Mode = FILES.Demo_Mode
        self.bohr2ang = 0.529177
        self.bohr2ang2 = self.bohr2ang * self.bohr2ang
        self.bohr2ang3 = self.bohr2ang2 * self.bohr2ang
        self.ang2bohr = 1/self.bohr2ang

    def InitBWlist(self, RawPar, logfilename):
        logfile = open(logfilename, "a")

        # if Use_AmGr_sel_crit is not present in input file, but some BWlist
        # is, Use_AmGr_sel_crit is set to True
        temp = 0
        usedict = {}
        SetNames = get_settings_names()
        for item in SetNames["settings_names_list_BWlist"]:
            try:
                temp += RawPar.input_parameters["Use_" + item]
            except Exception:
                usedict["Use_" + item] = False
            else:
                usedict["Use_" + item] = RawPar.input_parameters["Use_" + item]

        check = False
        if "Use_AmGroup_selection_criteria" not in RawPar.input_parameters:
            if not self.Demo_Mode:
                logfile.write(
                    "\nInput file does not specify whether to use Amide "
                    "group selection criteria. ")
            if temp != 0:
                self.Use_AmGroup_selection_criteria = True
                if not self.Demo_Mode:
                    logfile.write(
                        "As there are selection criteria present, these will "
                        "be used. If you don't want this, either "
                        "remove/comment out these selection criteria from the "
                        "input parameter file, or add the line "
                        "'Use_AmGroup_selection_criteria False' to the input "
                        "parameter file.")
            else:
                self.Use_AmGroup_selection_criteria = False
                if not self.Demo_Mode:
                    logfile.write(
                        "But as there are no selection criteria present, all "
                        "groups (without further selection) will be "
                        "calculated.")
            check = True
        elif RawPar.input_parameters["Use_AmGroup_selection_criteria"]:
            if temp == 0:
                logfile.write(
                    "\n\nAlthough selection criteria should be used, none are "
                    "specified. Therefore, they will not be used.")
                self.Use_AmGroup_selection_criteria = False
            else:
                self.Use_AmGroup_selection_criteria = True
                check = True
        else:
            logfile.write(
                "\n\nSelection criteria should not be used, so any specified "
                "criteria are disabled!")
            self.Use_AmGroup_selection_criteria = False
            self.Use_resname_blacklist = False
            self.Use_resname_whitelist = False
            self.Use_resnum_blacklist = False
            self.Use_resnum_whitelist = False

        if check:
            self.Use_resname_blacklist = usedict["Use_resname_blacklist"]
            self.Use_resname_whitelist = usedict["Use_resname_whitelist"]
            self.Use_resnum_blacklist = usedict["Use_resnum_blacklist"]
            self.Use_resnum_whitelist = usedict["Use_resnum_whitelist"]

        if self.Use_resname_blacklist:
            self.resname_blacklist = (
                RawPar.input_parameters["resname_blacklist"])
        if self.Use_resname_whitelist:
            self.resname_whitelist = (
                RawPar.input_parameters["resname_whitelist"])
        if self.Use_resnum_blacklist:
            self.resnum_blacklist = (
                RawPar.input_parameters["resnum_blacklist"])
        if self.Use_resnum_whitelist:
            self.resnum_whitelist = (
                RawPar.input_parameters["resnum_whitelist"])

        logfile.close()

    def InitFrames(self, RawPar, logfilename):
        # dealing with which frames to calculate
        temp = 0
        templist = [True, True, True]
        for num, item in enumerate(
                ["start_frame", "nFrames_to_calculate", "end_frame"]):
            if item not in RawPar.input_parameters:
                temp += 1
                templist[num] = False
        if temp == 0:
            if (
                RawPar.input_parameters["start_frame"]
                + RawPar.input_parameters["nFrames_to_calculate"]
                != RawPar.input_parameters["end_frame"]
            ):
                ErrorText = (
                    "\nContradictory parameters: start_frame, "
                    "nFrames_to_calculate and end_frame are all specified, "
                    "but contradict. Either specify just two, or make sure "
                    "they match. Quitting!")
                AIM_PC.warning(ErrorText, True, 0, logfilename, self)
            else:
                self.start_frame = RawPar.input_parameters["start_frame"]
                self.nFrames_to_calculate = (
                    RawPar.input_parameters["nFrames_to_calculate"])
                self.end_frame = RawPar.input_parameters["end_frame"]
        elif templist[0]:
            self.start_frame = RawPar.input_parameters["start_frame"]
            if templist[1]:
                self.nFrames_to_calculate = (
                    RawPar.input_parameters["nFrames_to_calculate"])
                self.end_frame = self.start_frame + self.nFrames_to_calculate
            elif templist[2]:
                self.end_frame = RawPar.input_parameters["end_frame"]
                self.nFrames_to_calculate = self.end_frame - self.start_frame
            else:
                self.end_frame = RawPar.default_parameters["end_frame"]
                self.nFrames_to_calculate = self.end_frame - self.start_frame
        elif templist[1]:
            self.nFrames_to_calculate = (
                RawPar.input_parameters["nFrames_to_calculate"])
            if templist[2]:
                self.end_frame = RawPar.input_parameters["end_frame"]
                self.start_frame = self.end_frame - self.nFrames_to_calculate
            else:
                self.start_frame = RawPar.default_parameters["start_frame"]
                self.end_frame = self.start_frame + self.nFrames_to_calculate
        elif templist[2]:
            self.end_frame = RawPar.input_parameters["end_frame"]
            self.start_frame = RawPar.default_parameters["start_frame"]
            self.nFrames_to_calculate = self.end_frame - self.start_frame
        else:
            self.start_frame = RawPar.default_parameters["start_frame"]
            self.end_frame = RawPar.default_parameters["end_frame"]
            self.nFrames_to_calculate = self.end_frame - self.start_frame

    def InitNSA(self, RawPar, logfilename):
        self.NSA_toggle = RawPar.input_parameters["NSA_toggle"]
        if self.NSA_toggle:
            for item in ["NSA_spheresize", "NSA_nframes"]:
                if item not in RawPar.input_parameters:
                    ErrorText = (
                        "\nInput file is missing the following required "
                        "parameter: " + item + "\nThis is only required as "
                        "NSA_toggle is set to true. Either specify this "
                        "parameter, or set NSA_toggle to false. Quitting!")
                    AIM_PC.warning(ErrorText, True, 0, logfilename, self)
            if RawPar.input_parameters["NSA_spheresize"] <= self.SphereSize:
                ErrorText = (
                    "\nThe spheresize for the NSA (NSA_spheresize) is too "
                    "small compared to the final sphere size for the "
                    "electrostatic interactions (SphereSize). NSA_spheresize "
                    "should be larger than SphereSize. Check the manual for "
                    "more information. Quitting!")
                AIM_PC.warning(ErrorText, True, 0, logfilename, self)
            self.NSA_spheresize = RawPar.input_parameters["NSA_spheresize"]
            self.NSA_nframes = RawPar.input_parameters["NSA_nframes"]
            if self.use_c_lib:
                self.NSA_spheresize_c = np.float32(self.NSA_spheresize)

    def FixContradictions(self, logfilename, FILES):
        # # not useful to run if no output is requested
        # if not self.output_bin and not self.output_txt:
        #     ErrorText = (
        #         "\nContradicting input parameters: according to these, "
        #         "neither binary nor text output files should be created.
        # With "
        #         "no output, there is no need to run! Quitting!")
        #     AIM_PC.warning(ErrorText, True, 0, logfilename, self)
        # if not self.output_Ham and not self.output_Dip and not
        # self.output_Pos:
        #     ErrorText = (
        #         "\nContradicting input parameters: according to these, "
        #         "neither hamiltonian nor dipole or position output files "
        #         "should be created. With no output, there is no need to
        # run! "
        #         "Quitting!")
        #     AIM_PC.warning(ErrorText, True, 0, logfilename, self)

        # as MDA has faulty COM's, it makes little/no sense to still do the
        # NSA optimization (use_c_lib does work)
        if self.replicate_orig_AIM:
            if self.NSA_toggle:
                ErrorText = (
                    "\nNSA was requested, but so was the replicating the "
                    "behaviour of the original AmideImaps.c script. However, "
                    "due to technical limitations of this replication, NSA "
                    "cannot be used, and will be disabled!")
                AIM_PC.warning(ErrorText, False, 1, logfilename, self)
                self.NSA_toggle = False

        # if the protein is not part of the selection, it is incorrect to
        # consider NN's!
        if "All" not in self.influencers and "Protein" not in self.influencers:
            if self.TreatNN:
                ErrorText = (
                    "\nTreatNN was requested, but as the influencers "
                    "selection does not contain any protein parts, the "
                    "influence of nearest neighbours should not be considered "
                    "either. Setting TreatNN to False!")
                AIM_PC.warning(ErrorText, False, 1, logfilename, self)
                self.TreatNN = False

        # the parameter atom_based_chainID is a must for amber files! (They
        # lack both molnums AND segnums...) So, if this parameter is not used,
        # see if any amber files are used.
        if (
            not self.atom_based_chainID
            and
            (
                any(FILES.topfile.upper().endswith(k) for k in [
                    "TOP", "PRMTOP", "PARM7"
                    ])
                or
                any(FILES.trjfile.upper().endswith(k) for k in [
                    "TRJ", "MDCRD", "CRDBOX", "NC", "NCDF"
                ])
            )
        ):
            ErrorText = (
                "\nThe usage of amber-style topology and/or trajectory files "
                "has been detected. These files require the protein chains "
                "to be identified using an atom-based method, which wasn't "
                "requested. Setting atom_based_chainID to True!"
            )
            AIM_PC.warning(ErrorText, False, 1, logfilename, self)
            self.atom_based_chainID = True

    def ValidateInput(self, logfilename):
        """
        Checks whether the choice for any parameter is valid - whether AIM
        understands it!
        """
        setdict = get_settings_allowed_values()
        problem = False

        for key in setdict:
            choice = getattr(self, key)

            if type(choice) == list:
                if any(
                    choicepart not in setdict[key] for choicepart in choice
                ):
                    problem = True
            else:
                if choice not in setdict[key]:
                    problem = True

            if problem:
                ErrorText = (
                    "Invalid choice for parameter " + key
                    + ": only the following choices are allowed:")
                for item in setdict[key]:
                    ErrorText = ErrorText + " " + str(item)
                ErrorText = ErrorText + "\nQuitting!"
                AIM_PC.warning(ErrorText, True, 0, logfilename, self)

    def ReadExtraMaps(self, FILES):
        # the map needs access to proper parameters.
        self.ExtraMaps = extra_map_reader(FILES, self)
        self.ExtraMaps_names = [v.name for v in self.ExtraMaps.values()]
        self.ExtraMaps_names.append("AmideBB")

        if self.replicate_orig_AIM:
            self.ExtraMaps_names = ["AmideSC", "AmideBB"]

        (
            self.CouplingMaps, self.Coupling_all_pairs
        ) = extra_coup_map_reader(FILES, self)

    def ValidateMaps(self, FILES):
        """
        Validates the parameters map_choice and Dipole_choice
        """

        av_Emaps = FILES.AllMaps.available_Emaps
        av_Dmaps = FILES.AllMaps.available_Dmaps + ["Torii"]

        if self.map_choice not in av_Emaps:
            ErrorText = (
                "Invalid choice of frequency map: only the following choices "
                "are allowed: "
                + " ".join(av_Emaps)
                + "\nQuitting!"
            )
            AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)

        if self.Dipole_choice not in av_Dmaps:
            ErrorText = (
                "Invalid choice of dipole map: only the following choices "
                "are allowed: "
                + " ".join(av_Dmaps)
                + "\nQuitting!"
            )
            AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)

    def ValidateExtraMaps(self, FILES):
        """
        Validates the choice of oscillators. This is a separate function, as
        this validation requires the extra maps to be loaded, while these
        extra maps in turn require that the RunPar class is finished.
        """

        for choice in self.oscillators:
            if choice not in self.ExtraMaps_names:
                ErrorText = (
                    "Invalid choice of oscillators: only the following choices"
                    " are allowed: "
                    + " ".join([i for i in self.ExtraMaps_names])
                    + "\nQuitting!"
                )
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)

    def finish_init(self, FILES):
        # max_time is requested in minutes, but used internally in seconds.
        self.max_time *= 60

        # Allows special selections of residues (as influencers)
        self.ResnamesReader(FILES)
        # specifies the names of important atoms in the different packages.
        self.AtnamesReader(FILES)

    def ResnamesReader(self, FILES):
        """
        Reads and parses the resnames file - what residue names are known, and
        what kind are they? These groups are mainly used for influencer
        selection in the input parameter file.
        """
        rawresnamesfile = open(FILES.resnamesfilefilename)
        all_resname_vars = []
        resnames = {}

        # Protected groups:
        resnames_Protein = [
            "ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN",
            "CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE",
            "TYR", "TRP"]
        resnames_protein_special = ["FOR", "ETA", "GL2"]

        if self.use_protein_specials:
            resnames_Protein += resnames_protein_special
        resnames["Protein"] = resnames_Protein

        for lineraw in rawresnamesfile:
            line = LineFormat(lineraw)

            # if the line is empty, skip!
            if len(line) == 0:
                continue

            if line[0] == 'add':
                if line[1] in ["Protein", "All"]:
                    ErrorText = (
                        "The resnames file contains a protected resname. "
                        "Please change it: " + line[1] + ". Quitting!")
                    AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)
                resnames[line[1]] = []
                while len(line) >= 3:
                    try:
                        resnames[line[1]].extend(resnames[line[2]])
                    except Exception:
                        ErrorText = (
                            "The definition of " + line[1] + " in the "
                            "resnamefile relies on " + line[2] + ", but this "
                            "resname has not been specified. Please make sure "
                            "to specify it on a line higher up than the line "
                            "where it is added to another group. Quitting!")
                        AIM_PC.warning(
                            ErrorText, True, 0, FILES.logfilename, self)
                    del line[2]
            elif line[0] in ["Protein", "All"]:
                ErrorText = ("The resnames file contains a protected resname. "
                             "Please change it: " + line[0] + ". Quitting!")
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)
            else:
                all_resname_vars.append(line[0])
                resnames[line[0]] = line[1:]

        all_resname_vars.append("Protein")

        if "Protein_ext" in all_resname_vars:
            for residue in resnames["Protein_ext"]:
                if residue not in resnames_Protein:
                    resnames["Protein"].append(residue)

        resnames["All"] = []
        for name in all_resname_vars:
            for residue in resnames[name]:
                if residue not in resnames["All"]:
                    resnames["All"].append(residue)

        influencers_resnames = []
        for seltype in self.influencers:
            if seltype not in resnames:
                ErrorText = (
                    "Invalid choice for the parameter 'influencers': a choice "
                    "may either be 'All' or 'Protein'. All other requested "
                    "groups should be specified in the resnames file "
                    + FILES.resnamesfilefilename + "\nQuitting!")
                AIM_PC.warning(ErrorText, True, 0, FILES.logfilename, self)
            for residue in resnames[seltype]:
                if residue not in influencers_resnames:
                    influencers_resnames.append(residue)

        self.ProtInSel = True
        for residue in resnames["Protein"]:
            if residue not in influencers_resnames:
                self.ProtInSel = False

        self.resnames_influencers = influencers_resnames
        self.resnames_Protein = resnames["Protein"]
        self.resnames_All = resnames["All"]

    def AtnamesReader(self, FILES):
        """
        Reads the atnames file. Depending on the MD package used, some atoms
        have different names. The atoms for which this is relevant are stored
        in the atnames file, along with the choice for them!
        """
        raw_atnamesfile = open(FILES.atnamesfilefilename)
        self.atnames = {}
        for lineraw in raw_atnamesfile:

            line = LineFormat(lineraw)

            # if the line is empty, skip!
            if len(line) < 2:
                continue

            self.atnames[line[0]] = line[1:]

    def SetOPLS(self, value):
        self.opls = value


class MapReader_old:  # deprecated
    """
    The class for reading built-in AIM maps (those in the sourcefiles
    directory).
    """
    def __init__(self, FILES, RunPar):
        self.EmapReader(FILES, RunPar)
        if RunPar.Dipole_choice == "Jansen":
            self.DmapReader(FILES, RunPar)
        else:
            self.Emap["Gen_Mu1_useG"] = False
            self.Emap["Pro_Mu1_useG"] = False
        self.NNmapReader(FILES.NN_Mapfilename)
        self.CouplingMapReader(
            FILES.TCC_Mapfilename, FILES.Tasumi_Mapfilename, RunPar)

        # self.OtherMapReader(FILES)

    def EmapReader(self, FILES, RunPar):
        function_warning = 0

        if RunPar.map_choice == "Tokmakoff":
            mapfilename = FILES.Tokmakoff_EMapfilename
        elif RunPar.map_choice == "Skinner":
            mapfilename = FILES.Skinner_EMapfilename
        elif RunPar.map_choice == "Jansen":
            mapfilename = FILES.Jansen_EMapfilename
        elif RunPar.map_choice == "Custom":
            mapfilename = FILES.Custom_EMapfilename

        EmapRaw, mapnames = self.EmapReader_unit(FILES, RunPar, mapfilename, 7)

        for mapname in ["Gen", "Pro", "SC"]:
            if mapname not in mapnames:
                ErrorText = (
                    "User opted for the " + RunPar.map_choice
                    + " map, saved as the file " + mapfilename
                    + "\nThe mapfile is missing a map of the name "
                    + mapname + ". Quitting!")
                AIM_PC.warning(
                    ErrorText, True, function_warning,
                    FILES.logfilename, RunPar)

        self.Emap = {}

        for mapname in mapnames:
            PEG = np.array(EmapRaw[mapname + "_PEG0"])
            PEGabs = abs(PEG)
            relevant_j = []
            for i in range(6):
                if np.amax(PEGabs[i, :]) > 0:
                    relevant_j.append(i)

            self.Emap[mapname+"_omega"] = EmapRaw[mapname+"_omega0"]
            self.Emap[mapname+"_P"] = PEG[:, 0]
            self.Emap[mapname+"_E"] = PEG[:, 1:4]
            self.Emap[mapname+"_G"] = PEG[:, 4:]

            self.Emap[mapname+"_relevant_j"] = relevant_j
            self.Emap[mapname+"_relevant_j_arr"] = np.array(
                relevant_j, dtype='int32')

            Gabs = abs(PEG[:, 4:])
            if np.amax(Gabs) > 0:
                self.Emap[mapname+"_useG"] = True
            else:
                self.Emap[mapname+"_useG"] = False

            if len(relevant_j) == 0:
                ErrorText = (
                    "User opted for the " + RunPar.map_choice
                    + " map, saved as the file " + mapfilename
                    + "\nThe mapfile contains nothing but zeros, which is a "
                    "problem. Please supply a map with more data. Quitting!")
                AIM_PC.warning(
                    ErrorText, True, function_warning,
                    FILES.logfilename, RunPar)

    def EmapReader_unit(self, FILES, RunPar, mapfilename, linespermap):
        function_warning = 0
        mapfile = open(mapfilename)

        linecounter = 0
        EmapRaw = {}
        mapname = ""
        mapnames = []
        for line in mapfile:
            line = LineFormat(line)
            if len(line) == 0:
                continue
            elif line[0] == "defmap":
                if len(line) == 1:
                    ErrorText = (
                        "User opted for the " + RunPar.map_choice
                        + " map, saved as the file " + mapfilename
                        + "\nThe keyword 'defmap' (indicating a new map) was "
                        "detected, but it was not followed by a name this map "
                        "should have! \nQuitting!")
                    AIM_PC.warning(
                        ErrorText, True, function_warning,
                        FILES.logfilename, RunPar)

                if linecounter != 0:
                    ErrorText = (
                        "User opted for the " + RunPar.map_choice
                        + " map, saved as the file " + mapfilename
                        + "\nA new map of the name " + line[1]
                        + " was declared, but the previous map is still "
                        "missing data! Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, function_warning,
                        FILES.logfilename, RunPar)

                linecounter = linespermap
                mapname = line[1]
            elif linecounter != 0 and linecounter % 7 == 0:
                linect = str((linecounter//7)-1)
                try:
                    omega = float(line[0])
                except Exception:
                    ErrorText = (
                        "User opted for the " + RunPar.map_choice
                        + " map, saved as the file " + mapfilename
                        + "\nThe value for omega for the " + line[1]
                        + " map was not recognised as a number.  Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, function_warning,
                        FILES.logfilename, RunPar)
                EmapRaw[mapname+"_omega"+linect] = omega
                EmapRaw[mapname+"_PEG"+linect] = []
                linecounter -= 1
            elif linecounter % 7 != 0:
                linect = str(linecounter//7)
                if len(line) != 10:
                    ErrorText = (
                        "User opted for the " + RunPar.map_choice
                        + " map, saved as the file " + mapfilename
                        + "\nThe values for P, E and G for the " + line[1]
                        + " map do not match the required format. Please "
                        "supply 10 numbers in total, separated by spaces. "
                        "Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, function_warning,
                        FILES.logfilename, RunPar)
                newline = []
                try:
                    for item in line:
                        newline.append(float(item))
                except Exception:
                    ErrorText = (
                        "User opted for the " + RunPar.map_choice
                        + " map, saved as the file " + mapfilename
                        + "\nThe values for P, E and G for the " + line[1]
                        + " aren't all numbers. Please specify numbers only. "
                        "Quitting!")
                    AIM_PC.warning(
                        ErrorText, True, function_warning,
                        FILES.logfilename, RunPar)
                EmapRaw[mapname+"_PEG"+linect].append(newline)

                if linecounter == 1:
                    mapnames.append(mapname)
                    mapname == ""

                linecounter -= 1

        mapfile.close()
        return EmapRaw, mapnames

    def DmapReader(self, FILES, RunPar):
        function_warning = 0

        mapfilename = FILES.Jansen_DMapfilename
        EmapRaw, mapnames = self.EmapReader_unit(
            FILES, RunPar, mapfilename, 21)

        for mapname in ["Gen_Mu1", "Pro_Mu1", "SC_Mu1"]:
            if mapname not in mapnames:
                ErrorText = (
                    "The dipole mapfile is missing a map of the name "
                    + mapname + ". Quitting!")
                AIM_PC.warning(
                    ErrorText, True, function_warning,
                    FILES.logfilename, RunPar)

        for mapname in mapnames:
            self.Emap[mapname+"_omega"] = []
            self.Emap[mapname+"_P"] = []
            self.Emap[mapname+"_E"] = []
            self.Emap[mapname+"_G"] = []
            useG = []
            for mudir in range(2, -1, -1):
                mudirs = str(mudir)
                PEG = np.array(EmapRaw[mapname+"_PEG"+mudirs])

                self.Emap[mapname+"_omega"].append(
                    EmapRaw[mapname+"_omega"+mudirs])
                self.Emap[mapname+"_P"].append(PEG[:, 0])
                self.Emap[mapname+"_E"].append(PEG[:, 1:4])
                self.Emap[mapname+"_G"].append(PEG[:, 4:])

                Gabs = abs(PEG[:, 4:])
                if np.amax(Gabs) > 0:
                    useG.append(True)
                else:
                    useG.append(False)

            if sum(useG) > 0:
                self.Emap[mapname+"_useG"] = True
            else:
                self.Emap[mapname+"_useG"] = False

        for mapname in ["Gen", "Pro", "SC"]:
            relevant_j = []
            for j in range(6):
                for P_E_G in ["P", "E", "G"]:
                    for direc in range(3):
                        name = mapname + "_Mu1_" + P_E_G
                        if P_E_G == "P":
                            maxval = np.amax(self.Emap[name][direc][j])
                        else:
                            maxval = np.amax(self.Emap[name][direc][j, :])
                        if maxval > 0 and j not in relevant_j:
                            relevant_j.append(j)

            # for j in self.Emap[mapname+"_relevant_j"]:
            #     if j not in relevant_j:
            #         relevant_j.append(j)
            relevant_j.sort()
            self.Emap[mapname+"_relevant_j_dipole"] = relevant_j
            self.Emap[mapname+"_Mu1_relevant_j_dipole_arr"] = np.array(
                relevant_j, dtype='int32')

    def NNmapReader(self, NN_Mapfilename):
        rawmapfile = open(NN_Mapfilename)
        self.NNmaps = {}
        for map_num in range(21):  # there are 21 maps - 3 for each situation!
            mapname = rawmapfile.readline().strip()
            for _ in range(13):
                _ = rawmapfile.readline()
            mapdata = np.genfromtxt(
                NN_Mapfilename, skip_header=map_num*14 + 1, max_rows=13)
            self.NNmaps[mapname] = mapdata

    def CouplingMapReader(self, TCC_Mapfilename, Tasumi_Mapfilename, RunPar):
        if (RunPar.coupling_choice == "TCC"
                or RunPar.NN_coupling_choice == "TCC"):
            rawmapfile = open(TCC_Mapfilename)
            temp = rawmapfile.readline().strip()
            self.TCC_4PiEps = float(temp)
            temp = rawmapfile.readline().strip().split()
            self.TCC_Gen_alpha = float(temp[0])
            self.TCC_Pro_alpha = float(temp[1])

            self.TCC_Gen_q = np.genfromtxt(
                TCC_Mapfilename, skip_header=2, max_rows=1)
            self.TCC_Pro_q = np.genfromtxt(
                TCC_Mapfilename, skip_header=3, max_rows=1)
            self.TCC_Gen_dq = np.genfromtxt(
                TCC_Mapfilename, skip_header=4, max_rows=1)
            self.TCC_Pro_dq = np.genfromtxt(
                TCC_Mapfilename, skip_header=5, max_rows=1)
            self.TCC_Gen_v = np.genfromtxt(
                TCC_Mapfilename, skip_header=7, max_rows=6)
            self.TCC_Pro_v = np.genfromtxt(
                TCC_Mapfilename, skip_header=14, max_rows=6)

            if RunPar.use_c_lib:
                # datatype management! save all data as c-friendly!
                self.TCC_4PiEps = np.float32(self.TCC_4PiEps)
                self.TCC_Gen_alpha = np.float32(self.TCC_Gen_alpha)
                self.TCC_Pro_alpha = np.float32(self.TCC_Pro_alpha)
                self.TCC_Gen_q = self.TCC_Gen_q.astype('float32')
                self.TCC_Gen_q_c = np.ctypeslib.as_ctypes(self.TCC_Gen_q)
                self.TCC_Pro_q = self.TCC_Pro_q.astype('float32')
                self.TCC_Pro_q_c = np.ctypeslib.as_ctypes(self.TCC_Pro_q)
                self.TCC_Gen_dq = self.TCC_Gen_dq.astype('float32')
                self.TCC_Gen_dq_c = np.ctypeslib.as_ctypes(self.TCC_Gen_dq)
                self.TCC_Pro_dq = self.TCC_Pro_dq.astype('float32')
                self.TCC_Pro_dq_c = np.ctypeslib.as_ctypes(self.TCC_Pro_dq)
                self.TCC_Gen_v_c = AIM_DC.ctype2d(self.TCC_Gen_v, 'float32')
                self.TCC_Pro_v_c = AIM_DC.ctype2d(self.TCC_Pro_v, 'float32')

        if RunPar.NN_coupling_choice == "Tasumi":
            self.Tasumi_Map = np.genfromtxt(Tasumi_Mapfilename, max_rows=13)


class ExtraMapReader:
    """
    The class for reading additional .map files supplied by the user. Every
    instance of this class corresponds to a single .map file. This makes it
    possible for users to store more information - they can save it as an
    attribute to that instance of this class.
    """
    # contains the contents of a single mapfile.
    def __init__(self, filename, FILES, RunPar):
        """
        This is called by the extra_map_reader function. It forms the base for
        the rest of this class - it calls every other function belonging to it.
        As long as functions are run successfully, the next one gets called,
        and certain attributes are assigned. Only when the .map file is read
        successfully, is the map actually added to the directory of extramaps.
        """
        self.filename = filename
        # self.extra_header = " for file " + os.path.basename(self.filename)
        self.extra_header = " of class ExtraMapReader"
        self.function_warning = 0
        rawfilecontents = self.opener(filename, FILES)
        if self.success:
            sortedcontents = self.grouper(rawfilecontents, FILES, filename)

        # as self.grouper can change the value of self.success, a new check
        # is necessary. (this happens more often)
        if self.success:
            self.reader(sortedcontents)

        if self.success:
            self.success = self.contentchecker(FILES, filename, RunPar)

        if self.success:
            self.use_G = self.functions["set_use_G"](FILES, RunPar)
            self.relevant_j = self.functions["set_relevant_j"](FILES, RunPar)
            self.relevant_j_arr = np.array(self.relevant_j, dtype='int32')
            self.references = self.functions["set_references"](
                self, FILES, RunPar)

            # if the inprot wasn't specified correctly in the ID section,
            # just set it to false.
            self.inprot = getattr(self, "inprot", False)

        AIM_PC.finwarning(self.function_warning, FILES.logfilename, FILES,
                          extra_header=self.extra_header)

    def opener(self, filename, FILES):
        """
        Extracts the contents from the .map file
        """
        try:
            file = open(filename)
            rawfilecontents = file.read()
        except Exception:
            ErrorText = (
                "\nFailed to open the file " + os.path.basename(filename) +
                " stored in the directory " + os.path.dirname(filename) +
                " . This problem occured while trying to read in the files "
                "containing the info on extramaps. This file will be skipped, "
                "but this could cause errors later. If you intend to use this "
                "file, please make sure it is present and formatted correctly."
            )
            self.function_warning = AIM_PC.warning(
                ErrorText, False, self.function_warning, FILES.logfilename,
                FILES, n=2, extra_header=self.extra_header)

            self.success = False
            return ""

        else:
            self.success = True
            return rawfilecontents

    def grouper(self, rawfilecontents, FILES, filename):
        """
        Reads the contents and looks for the headers. Saves the lines belonging
        to a header in a dictionary (with the header as key). A single value
        is a list of lines from the file. These lines have been parsed - all
        trailing whitespaces are removed, as well as empty lines and all text
        after '#'.
        """
        first = True
        sortedcontents = {}
        for line in rawfilecontents.split("\n"):

            # simple line formatting to ignore code after # (just like python)
            line = line.rstrip()
            line = line.split('#')[0]
            # hashpos = line.find("#")
            # if hashpos != -1:
            #     line = line[:hashpos]
            if len(line) == 0:
                continue

            linecheck = line.strip().strip("[] ")
            if linecheck in ["Identifiers", "Resname + Atnames", "Emap data",
                             "Dmap data", "Python Code", "References"]:
                if first:
                    section_name = linecheck
                    section_contents = []
                    first = False
                else:
                    sortedcontents[section_name] = section_contents
                    section_name = linecheck
                    section_contents = []

            else:
                section_contents.append(line)
        else:
            try:
                sortedcontents[section_name] = section_contents
            except Exception:
                ErrorText = (
                    "\nFailed to extract any information from the file "
                    + os.path.basename(filename) + " stored in the directory "
                    + os.path.dirname(filename) + ". This file will be "
                    "skipped, but this could cause errors later. If you "
                    "intend to use this file, please make sure that the "
                    "file is structured as specified in the manual."
                )

                self.function_warning = AIM_PC.warning(
                    ErrorText, False, self.function_warning, FILES.logfilename,
                    FILES, n=2, extra_header=self.extra_header)
                self.success = False
        return sortedcontents

    def reader(self, sortedcontents):
        """
        Takes the dictionary from grouper and reads/uses the contents. These
        are stored as attributes of the class. Each of the functions called
        here deals with that header in the .map file.
        """
        self.read_identifiers(sortedcontents)
        self.read_makeup(sortedcontents)
        self.read_Emapdata(sortedcontents)
        self.read_Dmapdata(sortedcontents)
        self.read_python_code(sortedcontents)
        self.read_references(sortedcontents)

    def read_identifiers(self, sortedcontents):
        datalist = sortedcontents.get("Identifiers", [])
        for line in datalist:
            line = line.split()
            if len(line) == 0:
                continue

            if line[0] == "name" and len(line) == 2:
                self.name = line[1]

            if line[0] == "ID" and len(line) == 2:
                try:
                    self.ID = int(line[1])
                except Exception:
                    pass

            if line[0] == "inprot" and len(line) == 2:
                if line[1] == "True":
                    self.inprot = True
                else:
                    self.inprot = False

    def read_makeup(self, sortedcontents):
        datalist = sortedcontents.get("Resname + Atnames", [])
        self.makeup = {}
        for line in datalist:
            line = line.split()

            if len(line) != 7:
                continue

            resnames = line[0].split(",")
            atnames = line[1:]

            for index, entry in enumerate(atnames):
                atnames[index] = entry.split(",")

            for resname in resnames:
                if resname in self.makeup:
                    self.makeup[resname].append(atnames)
                else:
                    self.makeup[resname] = [atnames]

    def read_python_code(self, sortedcontents):
        datalist = sortedcontents.get("Python Code", [])
        codestring = "\n".join(datalist)

        self.functions = {}
        # exec(codestring, None, self.functions)
        try:
            exec(codestring, None, self.functions)
        except Exception:
            self.exec_error = traceback.format_exc()
            # ErrorText = (
            #     "\n This occured while trying to read the code in the file "
            #     + self.filename)
            # print(dir(E))
            # raise type(E)(E.msg + ErrorText)

    def read_Emapdata(self, sortedcontents):
        datalist = sortedcontents.get("Emap data", [])
        temp_Emap = read_map_data(datalist)

        self.Emap = {}
        for mapname, data in temp_Emap.items():
            try:
                self.Emap[mapname] = np.asarray(data, dtype=np.float32)
            except Exception:
                pass

    def read_Dmapdata(self, sortedcontents):
        datalist = sortedcontents.get("Dmap data", [])
        temp_Dmap = read_map_data(datalist)

        self.Dmap = {}
        for mapname, data in temp_Dmap.items():
            try:
                self.Dmap[mapname] = np.asarray(data, dtype=np.float32)
            except Exception:
                pass

    def read_references(self, sortedcontents):
        datalist = sortedcontents.get("References", [])
        datastring = "\n".join(datalist)
        self.rawreferences = AIM_RM.readrefstring(datastring)

    def contentchecker(self, FILES, filename, RunPar):
        """
        Checks whether all required information is present. If anything
        substantial is missing, that instance is not saved for use. As all
        maps in the directory are read (but not necessarily used), this does
        not have to be an issue, so the user is warned, but AIM continues.

        It only is an issue when the user wants to use the map. In that case,
        AIM will still raise an error, as then, the input parameter file will
        ask for a non-existent map (that map failed somewhere here).
        """

        warnitems = []
        if not hasattr(self, "name"):
            warnitems.append("the unique name")

        if not hasattr(self, "ID"):
            warnitems.append("the ID number")

        if not hasattr(self, "makeup") or len(self.makeup) < 1:
            warnitems.append("residue name and/or atom names")

        if (
            not hasattr(self, "functions")
            or
            not all(k in self.functions for k in (
                "set_references",
                "set_use_G", "set_relevant_j", "local_finder", "pre_calc",
                "pre_frame", "calc_freq", "calc_dipole"
            ))
        ):
            warnitems.append("all required python functions")

        if "Ram" in RunPar.output_type:
            if "calc_raman" not in self.functions:
                warnitems.append("required raman function")

        gen_message_pre = "\nFailed to read "

        gen_message = (
            " for the file " +
            os.path.basename(filename) + " stored in the directory "
            + os.path.dirname(filename) + " . This file will be skipped, "
            "meaning the map cannot be used this calculation. If you intend "
            "to use this file, please make sure it is present and formatted "
            "correctly."
        )

        if len(warnitems) == 0:
            # no problems were found with the file!
            return True

        else:
            allwarns = ", ".join(warnitems)
            ErrorText = gen_message_pre + allwarns + gen_message
            if hasattr(self, "exec_error"):
                ErrorText += (
                    "\n\nReading the python code in this file failed raising "
                    "the following error:\n"
                    + self.exec_error
                )
            self.function_warning = AIM_PC.warning(
                ErrorText, False, self.function_warning, FILES.logfilename,
                FILES, n=2, extra_header=self.extra_header)
            return False


class CoupMapReader:
    """
    The class for reading additional .cmap files supplied by the user. Every
    instance of this class corresponds to a single .cmap file. This makes it
    possible for users to store more information - they can save it as an
    attribute to that instance of this class.
    """
    def __init__(self, filename, FILES, RunPar):
        """
        This is called by the extra_coup_map_reader function. It forms the base
        for the rest of this class - it calls every other function belonging to
        it. As long as functions are run successfully, the next one gets
        called, and certain attributes are assigned. Only when the .cmap file
        is read successfully, is the map actually added to the directory of
        coupmaps.
        """

        self.filename = filename
        self.extra_header = " of class CoupMapReader"
        self.function_warning = 0
        rawfilecontents = self.opener(filename, FILES)

        if self.success:
            sortedcontents = self.grouper(rawfilecontents, FILES, filename)

        if self.success:
            self.reader(sortedcontents)

        if self.success:
            self.success = self.contentchecker(FILES, filename)

        AIM_PC.finwarning(self.function_warning, FILES.logfilename, FILES,
                          extra_header=self.extra_header)

    def opener(self, filename, FILES):
        """
        Extracts the contents from the .map file
        """
        try:
            file = open(filename)
            rawfilecontents = file.read()
        except Exception:
            ErrorText = (
                "\nFailed to open the file " + os.path.basename(filename) +
                " stored in the directory " + os.path.dirname(filename) +
                " . This problem occured while trying to read in the files "
                "containing the info on coupling maps. This file will be "
                "skipped, but this could cause errors later. If you intend to "
                "use this file, please make sure it is present and formatted "
                " correctly."
            )
            self.function_warning = AIM_PC.warning(
                ErrorText, False, self.function_warning, FILES.logfilename,
                FILES, n=2, extra_header=self.extra_header)
            self.success = False
            return ""

        else:
            self.success = True
            return rawfilecontents

    def grouper(self, rawfilecontents, FILES, filename):
        """
        Reads the contents and looks for the headers. Saves the lines belonging
        to a header in a dictionary (with the header as key). A single value
        is a list of lines from the file. These lines have been parsed - all
        trailing whitespaces are removed, as well as empty lines and all text
        after '#'.
        """
        first = True
        sortedcontents = {}

        for line in rawfilecontents.split("\n"):

            # simple line formatting to ignore code after # (just like python)
            line = line.rstrip()
            hashpos = line.find("#")
            if hashpos != -1:
                line = line[:hashpos]

            linecheck = line.strip().strip("[] ")
            if linecheck in [
                "Identifiers", "Coupled OscIDs", "Cmap data", "Python Code",
                "References"
            ]:
                if first:
                    section_name = linecheck
                    section_contents = []
                    first = False
                else:
                    sortedcontents[section_name] = section_contents
                    section_name = linecheck
                    section_contents = []

            else:
                section_contents.append(line)
        else:
            try:
                sortedcontents[section_name] = section_contents
            except Exception:
                ErrorText = (
                    "\nFailed to extract any information from the file "
                    + os.path.basename(filename) + " stored in the directory "
                    + os.path.dirname(filename) + " . This file will be "
                    "skipped, but this could cause errors later. If you "
                    "intend to use this file, please make sure that the "
                    "file is structured as specified in the manual."
                )
                self.function_warning = AIM_PC.warning(
                    ErrorText, False, self.function_warning, FILES.logfilename,
                    FILES, n=2, extra_header=self.extra_header)
                self.success = False
        return sortedcontents

    def reader(self, sortedcontents):
        """
        Takes the dictionary from grouper and reads/uses the contents. These
        are stored as attributes of the class. Each of the functions called
        here deals with that header in the .map file.
        """
        self.read_identifiers(sortedcontents)
        self.read_coupled(sortedcontents)
        self.read_Cmapdata(sortedcontents)
        self.read_python_code(sortedcontents)
        self.read_references(sortedcontents)

    def read_identifiers(self, sortedcontents):
        datalist = sortedcontents.get("Identifiers", [])
        for line in datalist:
            line = line.split()
            if len(line) == 0:
                continue

            if line[0] == "ID" and len(line) == 2:
                try:
                    self.ID = int(line[1])
                except Exception:
                    pass

    def read_coupled(self, sortedcontents):
        datalist = sortedcontents.get("Coupled OscIDs", [])
        self.coupled = set()

        for line in datalist:
            line = line.split()

            if len(line) != 2:
                continue

            firstres = line[0].split(",")
            secres = line[1].split(",")

            try:
                firstres = [int(i) for i in firstres]
                secres = [int(i) for i in secres]
            except Exception:
                return

            for first in firstres:
                for sec in secres:
                    self.coupled.add((min(first, sec), max(first, sec)))

    def read_Cmapdata(self, sortedcontents):
        datalist = sortedcontents.get("Cmap data", [])
        temp_Cmap = read_map_data(datalist)

        self.Cmap = {}
        for mapname, data in temp_Cmap.items():
            try:
                self.Cmap[mapname] = np.asarray(data, dtype=np.float32)
            except Exception:
                pass

    def read_python_code(self, sortedcontents):
        datalist = sortedcontents.get("Python Code", [])
        codestring = "\n".join(datalist)

        self.functions = {}
        exec(codestring, None, self.functions)

    def read_references(self, sortedcontents):
        datalist = sortedcontents.get("References", [])
        datastring = "\n".join(datalist)
        self.references = AIM_RM.readrefstring(datastring)

    def contentchecker(self, FILES, filename):
        """
        Checks whether all required information is present. If anything
        substantial is missing, that instance is not saved for use. As all
        maps in the directory are read (but not necessarily used), this does
        not have to be an issue, so the user is warned, but AIM continues.

        The pair of groups that the failed map was supposed to treat will then
        be treated like any undefined pair of groups: no coupling, or just
        basic dipole-dipole coupling.
        """
        warnitems = []
        if not hasattr(self, "ID"):
            warnitems.append("the ID number")

        if not hasattr(self, "coupled") or len(self.coupled) < 1:
            warnitems.append("ID's of the groups to couple")

        if (
            not hasattr(self, "functions")
            or
            not all(k in self.functions for k in (
                "pre_calc", "prep_coupling", "calc_coupling"
            ))
        ):
            warnitems.append("all required python functions")

        gen_message_pre = "\nFailed to read "

        gen_message = (
            " for the file " +
            os.path.basename(filename) + " stored in the directory "
            + os.path.dirname(filename) + " . This file will be skipped, "
            "meaning the map cannot be used this calculation. If you intend "
            "to use this file, please make sure it is present and formatted "
            "correctly."
        )

        if len(warnitems) == 0:
            # no problems were found with the file!
            return True
        else:
            allwarns = ", ".join(warnitems)
            ErrorText = gen_message_pre + allwarns + gen_message
            self.function_warning = AIM_PC.warning(
                ErrorText, False, self.function_warning, FILES.logfilename,
                FILES, n=2, extra_header=self.extra_header)
            return False


def extra_map_reader(FILES, RunPar):
    """
    Called by RunParameters.__init__. Loops over the extramapdir, and tries to
    extract info from any file ending in .map. This is done with the
    ExtraMapReader class. When read successfully, those instances are then
    saved in a dict as values, with their ID as keys.
    """
    extra_mapdict = {}
    function_warning = 0
    # loop over all objects in extramaps directory
    for filename in os.listdir(FILES.extramapdir):
        # skip if the object is not a file
        if not os.path.isfile(os.path.join(FILES.extramapdir, filename)):
            continue
        # skip if the object is not a mapfile
        if not filename.endswith(".map"):
            continue
        extramap = ExtraMapReader(os.path.join(FILES.extramapdir, filename),
                                  FILES, RunPar)

        if not extramap.success:
            continue

        if (
            extramap.ID not in extra_mapdict
            and
            extramap.name not in [v.name for v in extra_mapdict.values()]
        ):
            extra_mapdict[extramap.ID] = extramap
        else:
            dupfilename = extra_mapdict[extramap.ID].filename
            ErrorText = (
                "\nThe map specified in the file " + os.path.basename(filename)
                + " located in the directory " + os.path.dirname(filename) +
                " has the same ID or name as the map specified in the file "
                + os.path.basename(dupfilename)
                + " located in the directory " + os.path.dirname(dupfilename) +
                " . Therefore, only the second of these files will be used, "
                "the first will be ignored. If you wish to use this file, "
                "please make sure it's name and ID are unique!"
            )
            function_warning = AIM_PC.warning(
                ErrorText, False, function_warning, FILES.logfilename,
                FILES)

    AIM_PC.finwarning(function_warning, FILES.logfilename, FILES)

    return extra_mapdict


def extra_coup_map_reader(FILES, RunPar):
    """
    Called by RunParameters.__init__. Loops over the extramapdir, and tries to
    extract info from any file ending in .cmap. This is done with the
    CoupMapReader class. When read successfully, those instances are then
    saved in a dict as values, with their ID as keys.
    """
    coup_mapdict = {}
    all_pairs = {}
    function_warning = 0

    # loop over all objects in extramaps directory
    for filename in os.listdir(FILES.extramapdir):
        # skip if the object is not a file
        if not os.path.isfile(os.path.join(FILES.extramapdir, filename)):
            continue
        # skip if the object is not a coupling file
        if not filename.endswith(".cmap"):
            continue
        coupmap = CoupMapReader(os.path.join(FILES.extramapdir, filename),
                                FILES, RunPar)

        if not coupmap.success:
            continue

        if coupmap.ID in coup_mapdict:
            dupfilename = coup_mapdict[coupmap.ID].filename
            ErrorText = (
                "\nThe map specified in the file " + os.path.basename(filename)
                + " located in the directory " + os.path.dirname(filename) +
                " has the same ID as the map specified in the file "
                + os.path.basename(dupfilename)
                + " located in the directory " + os.path.dirname(dupfilename) +
                " . Therefore, only the second of these files will be used, "
                "the first will be ignored. If you wish to use this file, "
                "please make sure it's ID is unique!"
            )
            function_warning = AIM_PC.warning(
                ErrorText, False, function_warning, FILES.logfilename, FILES)
            continue

        for pair in coupmap.coupled:
            if pair in all_pairs:
                dupfilename = coup_mapdict[all_pairs[pair]].filename
                ErrorText = (
                    "\nThe map specified in the file "
                    + os.path.basename(filename)
                    + " located in the directory "
                    + os.path.dirname(filename)
                    + " is coupling the groups " + str(pair[0])
                    + " and " + str(pair[1])
                    + ", just like the map specified in the file "
                    + os.path.basename(dupfilename)
                    + " located in the directory "
                    + os.path.dirname(dupfilename) +
                    ". As only type of coupling can be applied to any pair of "
                    "oscillators, only one coupling map can be defined for any"
                    " pair. Therefore, for this pair, only the map stored in "
                    + os.path.basename(dupfilename) +
                    " will be used, the other will be ignored. If you wish to "
                    "use the ignored file, please disable this pair in the "
                    "other file!"
                )
                function_warning = AIM_PC.warning(
                    ErrorText, False, function_warning, FILES.logfilename,
                    FILES)
            else:
                all_pairs[pair] = coupmap.ID

    AIM_PC.finwarning(function_warning, FILES.logfilename, FILES)

    return coup_mapdict, all_pairs


def read_map_data(datalist):
    """
    (support function for ExtraMapReader and CoupMapReader classes)
    Provided with a formatted data section, creates a dictionary with the data
    sorted by section.

    Input datalist: each line of data should be a separate entry in the list.
    The lines should already be formatted: no blank lines, no '#' and all
    text after '#'. This is because the input files of AIM (just like python)
    use '#' as an indication of non-computer text.

    output dict: if the correct data was provided by the user, there should be
    'defdata' keywords in the data. These keywords specify how that piece of
    data should be referred to. These names are used as keys in the output
    dict. The values are the arrays of numbers following the defdata keywords.
    """
    temp_map = {}
    mapdump = []
    dataname = ""
    for line in datalist:
        splitline = line.split()
        if len(splitline) == 0:
            continue
        if splitline[0] == "defdata":
            if dataname != "":
                temp_map[dataname] = mapdump
            mapdump = []
            if len(line) > 1:
                dataname = splitline[1]
            else:
                dataname = ""
        else:
            mapdump.append(splitline)
    else:
        if dataname != "":
            temp_map[dataname] = mapdump
    return temp_map


def get_settings_names():
    """
    This stores all valid parameter names that can be expected in the input and
    default parameter files. They are grouped by expected datatype of choices.
    """
    setnames = {}
    setnames["settings_names_list_float"] = [
        "SphereSize", "NSA_spheresize", "Scale_LRCoupling"]
    setnames["settings_names_list_int"] = [
        "Verbose", "Verbose_log", "start_frame", "nFrames_to_calculate",
        "end_frame", "max_time", "NSA_nframes"]
    setnames["settings_names_list_intBWlist"] = [
        "resnum_whitelist", "resnum_blacklist"]
    setnames["settings_names_list_strBWlist"] = [
        "resname_whitelist", "resname_blacklist"]
    setnames["settings_names_list_bool"] = [
        "use_c_lib", "profiler", "pngout", "replicate_orig_AIM", "NSA_toggle",
        "atom_based_chainID",
        "use_protein_specials", "Use_AmGroup_selection_criteria", "TreatNN"]
    setnames["settings_names_list_str"] = [
        "topfile", "trjfile", "sourcedir", "def_parfile", "referencefile",
        "resnamesfile", "atnamesfile", "libfile",
        "extramapdir", "outfilename", "outdipfilename",
        "outposfilename", "outparfilename", "outramfilename", "logfilename",
        "proffilename", "pngfilename", "apply_dd_coupling", "map_choice",
        "Dipole_choice",
        "coupling_choice", "NN_coupling_choice", "AtomPos_choice"]
    setnames["settings_names_list_str_list"] = [
        "output_format", "output_type", "influencers", "oscillators"
    ]
    setnames["settings_names_list_BWlist"] = (
        setnames["settings_names_list_intBWlist"]
        + setnames["settings_names_list_strBWlist"])

    setnames["settings_names_list_all"] = (
        setnames["settings_names_list_float"]
        + setnames["settings_names_list_int"]
        + setnames["settings_names_list_BWlist"]
        + setnames["settings_names_list_bool"]
        + setnames["settings_names_list_str"]
        + setnames["settings_names_list_str_list"]
    )

    return setnames


def get_settings_allowed_values():
    """
    For the parameters that only take certain choices, those allowed choices
    are stored here.
    """
    setdict = {}

    setdict["output_format"] = ["txt", "bin"]
    setdict["output_type"] = ["Ham", "Dip", "Pos", "Ram"]
    setdict["Verbose"] = [0, 1, 2, 3, 4]
    setdict["Verbose_log"] = [0, 1, 2, 3, 4]
    setdict["apply_dd_coupling"] = ["All", "None", "Same", "Different"]
    # setdict["map_choice"] = ["Tokmakoff", "Skinner", "Jansen", "Custom"]
    # setdict["Dipole_choice"] = ["Torii", "Jansen"]
    setdict["coupling_choice"] = ["TDCKrimm", "TDCTasumi", "TCC"]
    setdict["NN_coupling_choice"] = [
        "TDCKrimm", "TDCTasumi", "TCC", "Tasumi", "GLDP"]
    setdict["AtomPos_choice"] = ["C", "N", "O", "D"]

    return setdict


def read_parameter_file(parameter_filename, FILES, ftype, more_setname_str=[]):
    """
    The function that actually reads the input/default parameter files. Uses
    the get_settings_names() function to know what parameters to look for.
    Returns a dictionary with the parameter names and the formatted choices.
    Secondly, a dictionary with the raw input for each parameter is returned.
    The third object returned is a list with parameter names that were not
    recognised. This function does not verify whether choices are valid - the
    unrecognized parts are names of settings themselves.
    """

    SetNames = get_settings_names()
    function_warning = 0

    parameters_str = {}
    parameters_dict = {}
    unrecognised_parameters = []
    logfilename = FILES.logfilename

    try:
        parameter_file = open(parameter_filename)
    except FileNotFoundError:
        if ftype == "input":
            ErrorText = (
                "Could not open the file " + parameter_filename + " . This "
                "filename was specified in the command line while calling the "
                "program. Please make sure it is spelled correctly, and "
                "located in the correct directory. Quitting!"
            )
        else:
            # This is for when reading the default parameter file. That one
            # is already tested for before calling this function.
            ErrorText = "Some placeholder. You should never see this. FH.rpf"
        AIM_PC.warning(ErrorText, True, function_warning, logfilename, FILES)

    for lineraw in parameter_file:  # read parameter file line by line.

        line = LineFormat(lineraw)

        # if the line is empty, skip!
        if len(line) == 0:
            continue

        # accumulate all info stored on the line back into a string -> for
        # printing how the program reads the input file (doesn't print what
        # the program doesn't read)
        linestr = " ".join(line[1:])

        # save entered value as str
        parameters_str[line[0]] = linestr

        # now, interpret the information on the line to save as a parameter!

        # dealing with multiple strings
        if line[0] in SetNames["settings_names_list_str_list"]:
            parameters_dict[line[0]] = rpf_selection(
                line, parameter_filename, FILES, function_warning)

        # dealing with lines containing floats
        elif line[0] in SetNames["settings_names_list_float"]:
            out = rpf_float(line, parameter_filename, FILES, function_warning)
            if len(out) == 1:
                function_warning = out[0]
            else:
                function_warning = out[0]
                parameters_dict[line[0]] = out[1]

        # dealing with lines containing res black/whitelists
        elif line[0] in SetNames["settings_names_list_intBWlist"]:
            if len(line) == 1:
                parameters_dict["Use_" + line[0]] = False
            else:
                BWlistout = rpf_intBWlist(line, parameter_filename,
                                          FILES, function_warning)
                parameters_dict["Use_" + line[0]] = True
                parameters_dict[line[0]] = BWlistout

        elif line[0] in SetNames["settings_names_list_strBWlist"]:
            if len(line) == 1:
                parameters_dict["Use_" + line[0]] = False
            else:
                parameters_dict["Use_" + line[0]] = True
                parameters_dict[line[0]] = line[1:]

        # dealing with lines containing integers:
        elif line[0] in SetNames["settings_names_list_int"]:
            out = rpf_int(line, parameter_filename, FILES, function_warning)
            if len(out) == 1:
                function_warning = out[0]
            else:
                function_warning = out[0]
                parameters_dict[line[0]] = out[1]

        # dealing with lines containing booleans:
        elif line[0] in SetNames["settings_names_list_bool"]:
            out = rpf_bool(line, parameter_filename, FILES, function_warning)
            function_warning = out[0]
            parameters_dict[line[0]] = out[1]

        # dealing with lines containing strings:
        elif (line[0] in SetNames["settings_names_list_str"]
                or line[0] in more_setname_str):
            out = rpf_str(line, parameter_filename, FILES, function_warning)
            function_warning = out[0]
            parameters_dict[line[0]] = out[1]

        # if the specified parameter doesn't fit any of the known parameters,
        # print warning
        else:
            unrecognised_parameters.append(line[0])
            ErrorText = (
                "\nThe following specified parameter was not recognized: "
                + line[0] + "\nIt was found in the following file: "
                + parameter_filename + "\nDid you mean something else? The "
                "information will be ignored. Any missing parameter names "
                "will be appended from the default file. Check the run "
                "parameter file to see what the program runs with instead! "
                "You might want to stop this run, correct the input parameter "
                "file, and restart the calculation!")
            function_warning = AIM_PC.warning(
                ErrorText, False, function_warning, logfilename, FILES)

    parameter_file.close()
    AIM_PC.finwarning(function_warning, logfilename, FILES)

    return parameters_dict, parameters_str, unrecognised_parameters


def rpf_selection(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for the
    influencers and oscillators keywords. Returns a list of strings. Each item
    in the list is a group of choice.
    """
    if len(line) == 1:
        ErrorText = (
            "No choice detected for parameter " + line[0] + " in file "
            + parameter_filename
            + "\nThis entry will therefore be ignored, and the "
            "default setting will be used for this parameter instead. "
            "If you do not want this, cancel the calculation, and "
            "check the corresponding parameter file."
        )

        AIM_PC.warning(
            ErrorText, True, function_warning, FILES.logfilename, FILES, n=2)
    else:
        return line[1:]


def rpf_float(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for any
    parameter expecting a float. Returns the choice of the user as a float.
    """
    if len(line) == 1:
        ErrorText = (
            "No choice detected for parameter " + line[0] + " in file "
            + parameter_filename
            + "\nThis entry will therefore be ignored, and the "
            "default setting will be used for this parameter instead. "
            "If you do not want this, cancel the calculation, and "
            "check the corresponding parameter file.")
        function_warning = AIM_PC.warning(
            ErrorText, False, function_warning, FILES.logfilename, FILES, n=2)
        return function_warning
    try:
        returnval = float(line[1])
        return function_warning, returnval
    except Exception:
        ErrorText = (
            "Invalid choice detected for parameter " + line[0]
            + " in file " + parameter_filename
            + "\nPlease make sure that your choice is a number, with "
            "only a '.' to mark the decimal point, and no other "
            "special characters. Quitting!")
        AIM_PC.warning(
            ErrorText, True, function_warning, FILES.logfilename, FILES, n=2)


def rpf_intBWlist(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for the
    integer black-/whitelists. These accept dash notation, this is converted.
    Returns a list of integers.
    """
    BWlistin = line[1:]
    BWlistout = []

    for value in BWlistin:
        if "-" in value:
            cut = value.index("-")
            lowbound = value[:cut]
            upbound = value[cut+1:]

            try:
                lowbound = int(lowbound)
                upbound = int(upbound)
            except Exception:
                ErrorText = (
                    "Invalid choice detected for parameter " + line[0] +
                    " in file " + parameter_filename + "\nOne of the terms "
                    "was recognised as a range of residues (as it contained "
                    "a - character), but the numbers before and/or after it "
                    "were not recognised. Please make sure that the range is "
                    "defined by two numbers only. For example: 14-20 "
                    "\nQuitting!")
                AIM_PC.warning(
                    ErrorText, True, function_warning,
                    FILES.logfilename, FILES, n=2)

            if lowbound == upbound:
                BWlistout.append(lowbound)
            else:
                if lowbound > upbound:
                    temp = lowbound
                    lowbound = upbound
                    upbound = temp
                BWlistout.extend(range(lowbound, upbound+1, 1))

        else:  # no "-" in value
            try:
                BWlistout.append(int(value))
            except Exception:
                ErrorText = (
                    "Invalid choice detected for parameter "
                    + line[0] + " in file " + parameter_filename
                    + "\nPlease make sure that your choice is a "
                    "list of integer numbers (which contain "
                    "numbers only, no (special) characters), "
                    "and/or ranges (two integer numbers separated "
                    "by - instead of a space), separated by "
                    "spaces. Quitting!")
                AIM_PC.warning(
                    ErrorText, True, function_warning,
                    FILES.logfilename, FILES, n=2)

    return BWlistout


def rpf_int(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for any
    parameter expecting an integer. Returns the choice of the user as an
    integer.
    """
    if len(line) == 1:
        ErrorText = (
            "No choice detected for parameter " + line[0] + " in file "
            + parameter_filename + "\nThis entry will therefore be "
            "ignored, and the default setting will be used for this "
            "parameter instead. If you do not want this, cancel the "
            "calculation, and check the corresponding parameter file.")
        function_warning = AIM_PC.warning(
            ErrorText, False, function_warning, FILES.logfilename, FILES, n=2)
        return function_warning
    try:
        returnval = int(line[1])
        return function_warning, returnval
    except Exception:
        ErrorText = (
            "Invalid choice detected for parameter " + line[0]
            + " in file " + parameter_filename + "\nPlease make sure "
            "that your choice is a number, with no (special) "
            "characters. Quitting!")
        AIM_PC.warning(
            ErrorText, True, function_warning, FILES.logfilename, FILES, n=2)


def rpf_bool(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for any
    parameter expecting a boolean. Returns the choice of the user as a boolean.
    """
    if len(line) == 1:
        ErrorText = (
            "No choice detected for parameter " + line[0] + " in file "
            + parameter_filename + "\nThe choice for this entry will "
            "therefore be assumed to be 'True'. If you do not want "
            "this, cancel the calculation, and check the "
            "corresponding parameter file.")
        function_warning = AIM_PC.warning(
            ErrorText, False, function_warning, FILES.logfilename, FILES, n=2)
        line.append("True")
    if line[1] in ["True", "TRUE", "true"]:
        return function_warning, True
    elif line[1] in ["False", "FALSE", "false"]:
        return function_warning, False
    else:
        ErrorText = (
            "Invalid choice detected for parameter " + line[0]
            + " in file " + parameter_filename
            + "\nPlease make sure that your choice is either the "
            "text 'True', or 'False', without quotation marks. "
            "Quitting!")
        AIM_PC.warning(
            ErrorText, True, function_warning, FILES.logfilename, FILES, n=2)


def rpf_str(line, parameter_filename, FILES, function_warning):
    """
    Called by the read_parameter_file function. Parses the choice for any
    parameter expecting a string. Returns the choice of the user as a string.
    """
    if len(line) == 1:
        ErrorText = (
            "No choice detected for parameter " + line[0] + " in file "
            + parameter_filename + "\nThis entry will therefore be "
            "ignored, and the default setting will be used for this "
            "parameter instead. If you do not want this, cancel the "
            "calculation, and check the corresponding parameter file.")
        function_warning = AIM_PC.warning(
            ErrorText, False, function_warning, FILES.logfilename, FILES, n=2)
    returnval = line[1]
    return function_warning, returnval


def LineFormat(string):
    """
    General line formatting for when reading parameter files. Removes any
    leading/trailing whitespaces, and '#' along with everything following '#'.
    Then splits the line into a list, and returns the list.
    """
    # remove any white spaces/line breaks
    linestr = string.strip()

    # if there is a # somewhere in the line, ignore everything after it.
    linestr = linestr.split('#')[0]
    # check = -1
    # for index, value in enumerate(linestr):
    #     if value == "#":
    #         check = index
    #         break
    # if check != -1:
    #     linestr = linestr[:check]

    # split the line into a list
    line = linestr.split()

    return line


def NameEdit(string, FILES):
    """
    Output file names may contain the current date/time and AIM version. This
    is where the [now] and [version] syntax are replaced with their actual
    values.
    """
    if len(string) >= 5:
        for i in range(len(string)-4):
            if string[i:i+5] == "[now]":
                string = string[:i] + FILES.now_str + string[i+5:]
                break
    if len(string) >= 9:
        for i in range(len(string)-8):
            if string[i:i+9] == "[version]":
                string = string[:i] + FILES.script_version + string[i+9:]
    return string


def logwrite(logfile, pre, offset, post, logstrwidth=48):
    """
    Formats and writes the supplied data to the log file. All is meant for a
    single line in the logfile.

    pre: the text for the first column. Starts after 'offset' spaces.
    offset: how many spaces should come before the pre text.
    post: the text that should be in the second column.
    logstrwidth: How many characters the first column should take up. Or, after
    how many characters the second column (and thus the post text) should
    start. Default = 48.
    """
    logfile.write("\n" + " "*offset)
    logfile.write(format(pre, "<" + str(logstrwidth - offset)) + post)


def parwrite(parfile, pre, offset, post, parstrwidth=48):
    """
    Formats and writes the supplied data to the parameter file. All is meant
    for a single line in the parameter file.

    pre: the text for the first column. Starts after 'offset' spaces.
    offset: how many spaces should come before the pre text.
    post: the text that should be in the second column.
    logstrwidth: How many characters the first column should take up. Or, after
    how many characters the second column (and thus the post text) should
    start. Default = 48.

    The function is very similar to logwrite, but deals differently with
    parameters starting with 'raw'. This is the case for specified file
    locations. The user specifies them relative to a certain directory, or
    in an absolute manner. In general, relative is easier to understand.
    However, for AIM, internally, absolute is a must. Therefore, the user's
    choice is saved as 'raw', so it can be printed as the user gave it. The
    'raw'-less version becomes the absolute, internal one.
    """
    parfile.write("\n" + " "*offset)
    if pre[-3:] == "raw":
        pre = pre[:-3]
    parfile.write(format(pre, "<" + str(parstrwidth - offset)) + str(post))


def bulkparwrite(parfile, offset, itemslist, ParStore, parstrwidth=48):
    """
    For writing many parameters to the parameter file at once. Makes use of
    the parwrite function.
    """
    for item in itemslist:
        try:
            if type(getattr(ParStore, item)) is list:
                parwrite(parfile, item, offset, "")
                for ititem in getattr(ParStore, item):
                    parfile.write(str(ititem) + " ")
            else:
                parwrite(parfile, item, offset, str(getattr(ParStore, item)))
        except Exception:
            pass


def dummy_func(RunPar):
    """
    this is just to use the otherwise unused imports at the top, to avoid
    the warnings raised during editing this script.
    these are needed as an import, as the .map and the .cmap files might
    require them.
    """

    vec = AIM_SO.threesort_KvA([3, 4, 2])
    vec = np.array(vec)
    vec /= AIM_MF.vec3_len(vec)
    vec = AIM_PF.Dipole_Torii(vec, vec)

    met1, met2 = AIM_CF.DetCalcMeth(RunPar)
    metarr = [met1]
    metarr.append(met2)

    timer = AIM_TK.Timer()
    timer.Endtime()
