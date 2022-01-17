
def pytest():
    global success
    try:
        import sys
    except Exception:
        print(
            "Cannot import the sys module. This should be included during the "
            "installation of python itself. Please make sure the installation "
            "has been completed correctly, as python (and thus AIM) cannot "
            "function like this. "
        )
        quit()
    else:
        print("Succesfully imported the Sys module for python version check")

    try:
        pyvers = sys.version_info
    except Exception:
        print(
            "Could not extract the python version using sys. Please make sure "
            "you are using at least python version 3.8.3!"
        )
        quit()
    else:
        print(
            "Succesfully retrieved python version number using the Sys module"
        )

    if pyvers[0] != 3:
        print(
            "Unsuitable version of python installation detected! Please make "
            "sure to run this script (and all others within AIM) using at "
            "least python version 3.8.3!"
        )
        quit()
    elif pyvers[1] < 8:
        print(
            "An old version of python has been detected! AIM has been tested "
            "using python versions 3.8.3 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False
    elif pyvers[1] == 8 and pyvers[2] < 3:
        print(
            "An old version of python has been detected! AIM has been tested "
            "using python versions 3.8.3 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False

    print(
        "Detected python version: " + str(pyvers[0]) + "." + str(pyvers[1])
        + "." + str(pyvers[2])
    )

    print("Completed test of python version.\n")


def numpytest():
    global success
    try:
        import numpy as np
    except Exception:
        print(
            "Could not find the numpy module. Make sure it is installed "
            "correctly. Instructions for this can be found on the numpy "
            "web pages."
        )
        quit()
    else:
        print("Numpy has been imported correctly!")

    try:
        npvers = [int(i) for i in np.__version__.split(".")]
    except Exception:
        print(
            "Could not retrieve version number of numpy installation. Please "
            "make sure numpy is installed correctly. Instructions for this "
            "can be found on the numpy web pages."
        )
        quit()
    else:
        print("Successfully retrieved numpy version number")

    if npvers[0] != 1:
        print(
            "An old version of numpy has been detected! AIM has been tested "
            "using numpy versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False
    elif npvers[1] < 18:
        print(
            "An old version of numpy has been detected! AIM has been tested "
            "using numpy versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False
    elif npvers[1] == 18 and npvers[2] < 5:
        print(
            "An old version of numpy has been detected! AIM has been tested "
            "using numpy versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False

    print(
        "Detected numpy version: " + str(npvers[0]) + "." + str(npvers[1])
        + "." + str(npvers[2])
    )

    print("Completed test of numpy version.\n")


def numbatest():
    global success
    try:
        import numba as nb
    except Exception:
        print(
            "Could not find the numba module. Make sure it is installed "
            "correctly. Instructions for this can be found on the numba "
            "web pages."
        )
        quit()
    else:
        print("Numba has been imported correctly!")

    try:
        nbvers = [int(i) for i in nb.__version__.split(".")]
    except Exception:
        print(
            "Could not retrieve version number of numba installation. Please "
            "make sure numba is installed correctly. Instructions for this "
            "can be found on the numba web pages."
        )
        quit()
    else:
        print("Successfully retrieved numba version number")

    if nbvers[0] != 0:
        print(
            "An old version of numba has been detected! AIM has been tested "
            "using numba versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False
    elif nbvers[1] < 50:
        print(
            "An old version of numba has been detected! AIM has been tested "
            "using numba versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False
    elif nbvers[1] == 50 and nbvers[2] < 1:
        print(
            "An old version of numba has been detected! AIM has been tested "
            "using numba versions 1.18.5 or newer. It is recommended that you "
            "use a newer version, proceed at your own risk!"
        )
        success = False

    print(
        "Detected numba version: " + str(nbvers[0]) + "." + str(nbvers[1])
        + "." + str(nbvers[2])
    )

    print("Completed test of numba version.\n")


def MDAtest():
    global success
    try:
        import MDAnalysis as MDA
    except Exception:
        print(
            "Could not find the MDAnalysis module. Make sure it is installed "
            "correctly. Instructions for this can be found on the MDAnalysis "
            "web pages."
        )
        quit()
    else:
        print("MDAnalysis has been imported correctly!")

    try:
        mdavers = [int(i) for i in MDA.__version__.split(".")]
    except Exception:
        print(
            "Could not retrieve version number of MDAnalysis installation. "
            "Please make sure MDAnalysis is installed correctly. Instructions "
            "for this can be found on the MDAnalysis web pages."
        )
        quit()
    else:
        print("Successfully retrieved MDAnalysis version number")

    if mdavers[0] != 2:
        print(
            "An old version of MDAnalysis has been detected! AIM has been "
            "tested using both MDAnalysis versions 1.0.0 and 2.0.0. While "
            "AIM could work using 1.0.0, it is recommended that you "
            "use a newer version, as the old version can only read in files "
            "created using older versions of MD software. "
            "Proceed at your own risk!"
        )
        success = False

    print(
        "Detected MDAnalysis version: " + str(mdavers[0]) + "."
        + str(mdavers[1]) + "." + str(mdavers[2])
    )

    print("Completed test of MDAnalysis version.\n")


print("Starting test of python and module versions!\n")

success = True

pytest()
numpytest()
numbatest()
MDAtest()

if success:
    print(
        "Python itself and all modules required by AIM were found to be "
        "present and of the correct version. You can now go ahead and use "
        "AIM. Have fun!"
    )
else:
    print(
        "While python itself and all modules required by AIM were found to be "
        "present, one or more of them are an outdated version. Please "
        "consider updating these. As AIM has not been tested in this "
        "configuration, proceed at your own risk!"
    )
