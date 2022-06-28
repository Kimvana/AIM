# standard lib imports
import os
import sys
import shutil

# my lib imports
# import .sourcecode.FileHandler as AIM_FH
# import .sourcecode.TimeKeeping as AIM_TK

from .sourcecode import FileHandler as AIM_FH
from .sourcecode import TimeKeeping as AIM_TK

def run(callcommand):
    TIMER = AIM_TK.Timer()

    # very basic: finds hard-coded log, output and sourcefiles dirs. Cleans up
    # the
    # log and output dirs. Determines the start time and OS used for this run.
    # FILES = AIM_FH.FileLocations(sys.path[0], "V9-11", "9.11")
    FILES = AIM_FH.FileLocations(__file__, "V1-0-1", "1.0.1")

    # A lot of parameter processing. Finds the requested and default
    # parameters,
    # solves any issues/clashes with the chosen parameter settings. Checks if
    # all
    # requested files are present. Writes choices to log file and creates the
    # parameter output file.
    print(sys.argv)
    print(type(sys.argv))
    RunPar = FILES.GenerateRunParameters(callcommand)


    # Prints the first logo to terminal, prints version, starts profiler. The
    # returned WS contains the topology/trajectory as read from the file, with
    # minor changes - just those to make WS independent of the MD package used.
    WS = FILES.InitProgram(RunPar)

    # Actually prepares the system;
    WS.InitSys(FILES, RunPar)
    WS.RunCalc(TIMER, FILES, RunPar)

    WS.FinishCalc(TIMER, FILES, RunPar)


def setup(target_folder):
    AIMdir = os.path.abspath(os.path.dirname(__file__))
    shutil.copytree(
        os.path.join(AIMdir, 'extra_maps'),
        os.path.join(target_folder, 'extra_maps')
    )
    shutil.copytree(
        os.path.join(AIMdir, 'sourcefiles'),
        os.path.join(target_folder, 'sourcefiles')
    )


def main():
    callcommand = sys.argv
    if callcommand[1] == "run":
        run(callcommand[1:])
    elif callcommand[1] == "setup":
        setup(callcommand[2])


if __name__ == "__main__":
    main()
