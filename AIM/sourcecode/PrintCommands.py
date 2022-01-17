import sys


def vprint(threshold, texttoprint, logfilename, RunPar):
    """
    AIM's alternative to the python print function. Prints both to the console
    and the log file, but only if the verbose setting reaches the threshold.
    """

    if getattr(RunPar, "Demo_Mode", False):
        # Start every line with the text "Demo Mode: ", to make it very clear
        # the program is using the included data, not the user's own data.
        for i in range(len(texttoprint)-1, -1, -1):
            if texttoprint[i] == "\n":
                texttoprint = (texttoprint[:i+1] + "Demo Mode: "
                               + texttoprint[i+1:])
        if len(texttoprint) < 9 or texttoprint[:9] != "Demo Mode":
            texttoprint = "Demo Mode: " + texttoprint

    Verbose_log = getattr(RunPar, "Verbose_log", 4)
    if Verbose_log >= threshold:
        try:
            logfile = open(logfilename, "a")
        except Exception:
            logfile = open(logfilename, "w")
        logfile.write(texttoprint + "\n")
        logfile.close()

    Verbose = getattr(RunPar, "Verbose", 4)
    if Verbose >= threshold:
        print(texttoprint)


def vprintl(threshold, listtoprint, logfilename, RunPar):
    """
    AIM's alternative to the python print function when multiple arguments
    are used (kind of like a list)
    """
    texttoprint = ""
    for item in listtoprint:
        texttoprint += str(item) + " "
    vprint(threshold, texttoprint, logfilename, RunPar)


def warning(
    message_text, forcestop, function_warning, logfilename, RunPar, n=1,
    extra_header=""
):
    """
    Prints the warning messages. If that function didn't give a warning before,
    start the warning print with a header containing the function name.
    """
    headerin = "\n\n====================\n"
    headerout = "\n===================="

    frame = sys._getframe(n)
    test = frame.f_code.co_name

    if test == "<module>":
        test = "AIM.py"
    # elif test == "__init__":
    #     print("__init__ found!")
    #     exinf = sys.exc_info()
    #     print(exinf)
    #     print(exinf[1].__traceback__)
    #     tb = exinf[2]
    #     print(tb)

    #     while True:
    #         print("---", tb)
    #         for i in dir(tb):
    #             print(i, ':', getattr(tb, i))
    #         tb = tb.tb_next
    #         if not tb:
    #             break

    #     quit()
        # func = types.FunctionType(frame.f_code, frame.f_globals)
        # print(func)
        # print(frame)
        # print(frame.f_code.__class__)
        # print(func.__self__.__class__.__name__)
    header = headerin + test + extra_header + headerout

    if message_text[0] != "\n":
        message_text = "\n" + message_text

    BL_presc = hasattr(RunPar, "backlog")

    if forcestop:
        if BL_presc:
            for threshold, message in RunPar.backlog:
                vprint(threshold, message, logfilename, RunPar)
        warn_core(
            0, header, message_text, logfilename, RunPar, function_warning
        )
        quit()

    if BL_presc:
        if function_warning == 0:
            RunPar.backlog.append([0, header])
            function_warning = 1
        RunPar.backlog.append([0, message_text])
    else:
        function_warning = warn_core(
            0, header, message_text, logfilename, RunPar, function_warning
        )

    return function_warning


def warn_core(threshold, header, message_text, logfilename, RunPar,
              function_warning):
    """
    Used by the function 'warning', to avoid copying too much code.
    """
    if function_warning == 0:
        vprint(threshold, header, logfilename, RunPar)
        function_warning = 1
    vprint(threshold, message_text, logfilename, RunPar)
    return function_warning


def finwarning(function_warning, logfilename, RunPar, n=1, extra_header=""):
    """
    The complement to the warning function. Like the warning function prints
    a header before the first error, this can print a footer after the last.
    """
    if function_warning == 1:
        headerin = "\n====================\n"
        headerout = "\n===================="
        test = sys._getframe(n).f_back.f_code.co_name
        if test == "<module>":
            test = "AIM.py"
        header = headerin + "End of " + test + extra_header + headerout

        if hasattr(RunPar, "backlog"):
            RunPar.backlog.append([0, header])
        else:
            vprint(0, header, logfilename, RunPar)
