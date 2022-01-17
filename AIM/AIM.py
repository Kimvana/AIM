# standard lib imports
import sys

# my lib imports
import sourcecode.FileHandler as AIM_FH
import sourcecode.TimeKeeping as AIM_TK


TIMER = AIM_TK.Timer()

# very basic: finds hard-coded log, output and sourcefiles dirs. Cleans up the
# log and output dirs. Determines the start time and OS used for this run.
FILES = AIM_FH.FileLocations(sys.path[0], "V1-0-0", "1.0.0")

# A lot of parameter processing. Finds the requested and default parameters,
# solves any issues/clashes with the chosen parameter settings. Checks if all
# requested files are present. Writes choices to log file and creates the
# parameter output file.
RunPar = FILES.GenerateRunParameters(sys.argv)


# Prints the first logo to terminal, prints version, starts profiler. The
# returned WS contains the topology/trajectory as read from the file, with
# minor changes - just those to make WS independent of the MD package used.
WS = FILES.InitProgram(RunPar)

# Actually prepares the system;
WS.InitSys(FILES, RunPar)
WS.RunCalc(TIMER, FILES, RunPar)

WS.FinishCalc(TIMER, FILES, RunPar)
