# standard lib imports
import sys

# my lib imports
from . import DataConversion as AIM_DC


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


def finprint_composition(FILES, RunPar, WS):
    """
    At the end of the calculation, a lot of information is printed to the
    console. This function manages the prints about the system: how many groups
    of what type were found in the system?
    """

    nchains = len(WS.fincalcprints["molnum_prot_chains"])
    vprintl(1, [
        "amount of chains: ", nchains,
        "\n"], FILES.logfilename, RunPar)
    nBBrawall = 0
    nBBall = 0
    counter = 0
    totaloftype = {}
    accounted_groups = set()
    for molnum in WS.fincalcprints["molnum_prot_chains"]:
        if "AmideBB" in RunPar.oscillators:
            nBBraw = WS.fincalcprints["nBBraw" + str(molnum)]
            nBB = WS.fincalcprints["nBBused" + str(molnum)]
            tothere = nBB
        else:
            tothere = 0
        vprint(
            2, "Chain " + str(counter) + ": ", FILES.logfilename, RunPar)
        if "AmideBB" in RunPar.oscillators:
            vprint(
                3, "Amide groups in backbone: " + str(nBBraw),
                FILES.logfilename, RunPar)
            vprint(
                2, "Amide groups in backbone considered: " + str(nBB),
                FILES.logfilename, RunPar)
        for ID in WS.fincalcprints["OscID_present"]:
            if not RunPar.ExtraMaps[ID].inprot:
                continue
            nOsc = 0
            for OscGroup in range(WS.res_desired_len):
                if (
                    WS.OscID[OscGroup] == ID
                    and
                    WS.molnums[WS.AllOscGroups[OscGroup, 6]] == molnum
                ):
                    nOsc += 1
                    accounted_groups.add(OscGroup)
            if ID in totaloftype:
                totaloftype[ID] += nOsc
            else:
                totaloftype[ID] = nOsc
            tothere += nOsc
            IDname = RunPar.ExtraMaps[ID].name
            vprint(
                2, "Groups of type " + IDname + " considered: "
                + str(nOsc), FILES.logfilename, RunPar)
        vprint(
            2, "Total amount of groups considered: "
            + str(tothere) + "\n", FILES.logfilename, RunPar)
        counter += 1
        if "AmideBB" in RunPar.oscillators:
            nBBrawall += nBBraw
            nBBall += nBB

    # after looking at protein chains, now consider the groups not part of
    # a protein chain
    vprint(2, "Groups outside of protein: ", FILES.logfilename, RunPar)
    tothere = 0
    for ID in WS.fincalcprints["OscID_present"]:
        nOsc = 0
        for OscGroup in range(WS.res_desired_len):
            if (
                WS.OscID[OscGroup] == ID
                and
                OscGroup not in accounted_groups
            ):
                nOsc += 1
        if ID in totaloftype:
            totaloftype[ID] += nOsc
        else:
            totaloftype[ID] = nOsc
        tothere += nOsc
        IDname = RunPar.ExtraMaps[ID].name
        vprint(
            2, "Groups of type " + IDname + " considered: "
            + str(nOsc), FILES.logfilename, RunPar)
    vprint(
        2, "Total amount of groups considered: "
        + str(tothere) + "\n", FILES.logfilename, RunPar)
    vprint(1, "\nSumming up:", FILES.logfilename, RunPar)
    if "AmideBB" in RunPar.oscillators:
        vprint(
            2, "Amide groups in backbone: " + str(nBBrawall),
            FILES.logfilename, RunPar)
        vprint(
            1, "Amide groups in backbone considered: " + str(nBBall),
            FILES.logfilename, RunPar)
    for ID in WS.fincalcprints["OscID_present"]:
        nOsc = totaloftype[ID]
        IDname = RunPar.ExtraMaps[ID].name
        vprint(
            2, "Groups of type " + IDname + " considered: "
            + str(nOsc), FILES.logfilename, RunPar)
    vprint(
        1, "Total amount of considered groups: "
        + str(WS.res_desired_len), FILES.logfilename, RunPar)

    if "AmideBB" in RunPar.oscillators and nBBall > 0:
        totaloftype[0] = nBBall
    return totaloftype


def finprint_frames(FILES, RunPar, WS):
    """
    At the end of the calculation, a lot of information is printed to the
    console. This function manages the prints about the treated frames. How
    many were available, and how many were actually calculated?
    """

    vprint(1, "\n\nCalculated frames:", FILES.logfilename, RunPar)
    vprint(
        1, "\nAnalyzed frames: " + str(RunPar.start_frame) + " to " +
        str(RunPar.end_frame-1), FILES.logfilename, RunPar)
    vprint(
        1, "Available frames: 0 to " + str(len(WS.trj)-1),
        FILES.logfilename, RunPar)


def finprint_time(FILES, RunPar, WS, TIMER):
    """
    At the end of the calculation, a lot of information is printed to the
    console. This function manages the prints about the calculation time:
    how long did the calculation take?
    """

    vprint(1, "\n\nCalculation time:", FILES.logfilename, RunPar)
    TIMER.Endtime()
    vprintl(1, [
        "\ninitialization:", round(TIMER.afterinit0, 3),
        "s\ntime per frame",
        round((TIMER.totheavy)/(
            WS.lastCalcFrame - RunPar.start_frame + 1), 3),
        "s\ntotal:", round((TIMER.endtime0), 3), "s =",
        AIM_DC.timestring(round((TIMER.endtime0), 3))
        ], FILES.logfilename, RunPar)


def finprint_references(FILES, RunPar, WS, totaloftype):
    """
    At the end of the calculation, a lot of information is printed to the
    console. This function manages the prints about the references: which
    papers should the user reference when they use this calculation for their
    own work?
    """

    vprint(1, "\n\nRelevant references:", FILES.logfilename, RunPar)

    # dict of strings of refs + lists of things using them!
    allrefs = {}
    AIMrefs = []

    # collect all ID's on which a certain calculation is performed
    freqID = []
    dipID = []
    ramID = []
    coupID = []

    if "Dip" in RunPar.output_type:
        for key in totaloftype.keys():
            if key not in dipID:
                dipID.append(key)

    if "Ham" in RunPar.output_type:
        coupID = WS.coup_used.keys()
        for key in totaloftype.keys():
            if key not in freqID:
                freqID.append(key)

        # for TDCgen, the calculated dipole moments are used.
        if 1 in WS.coup_used:
            for OscGroup in WS.coup_used[1]:
                if WS.OscID[OscGroup] not in dipID:
                    dipID.append(WS.OscID[OscGroup])

    if "Ram" in RunPar.output_type:
        for key in totaloftype.keys():
            if key not in ramID:
                ramID.append(key)

    if 0 in freqID or 0 in dipID:
        usedmaptypes = ["mapGen"]
        found = False
        for OscGroup in range(WS.res_desired_len):
            if WS.PrePro[OscGroup] == 1:
                found = True
                break
        if found:
            usedmaptypes.append("mapPro")

    for reference in FILES.OtherRefs:
        if "AIMcode" in reference.AIMnotes:
            refstr = str(reference)
            AIMrefs.append(refstr)

    for ID in freqID:
        if ID == 0:
            for reference in WS.AllMaps.Emap.references:
                if any("E"+x in reference.AIMnotes for x in usedmaptypes):
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][0].append("AmideBB")
                    else:
                        allrefs[refstr] = [
                            ["AmideBB"],
                            [],
                            [],
                            []
                        ]
        else:
            for entry in RunPar.ExtraMaps[ID].references:
                if "frequency" in entry[1]:
                    refstr = str(entry[0])
                    if refstr in allrefs:
                        allrefs[refstr][0].append(RunPar.ExtraMaps[ID].name)
                    else:
                        allrefs[refstr] = [
                            [RunPar.ExtraMaps[ID].name],
                            [],
                            [],
                            []
                        ]

    for ID in dipID:
        if ID == 0:
            if RunPar.Dipole_choice == "Torii":
                for reference in FILES.OtherRefs:
                    if "DmapTorii" in reference.AIMnotes:
                        refstr = str(reference)
                        if refstr in allrefs:
                            allrefs[refstr][1].append("AmideBB")
                        else:
                            allrefs[refstr] = [
                                [],
                                ["AmideBB"],
                                [],
                                []
                            ]
                continue
            for reference in WS.AllMaps.Dmap.references:
                if any("D"+x in reference.AIMnotes for x in usedmaptypes):
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][1].append("AmideBB")
                    else:
                        allrefs[refstr] = [
                            [],
                            ["AmideBB"],
                            [],
                            []
                        ]
        else:
            for entry in RunPar.ExtraMaps[ID].references:
                if "dipole" in entry[1]:
                    refstr = str(entry[0])
                    if refstr in allrefs:
                        allrefs[refstr][1].append(RunPar.ExtraMaps[ID].name)
                    else:
                        allrefs[refstr] = [
                            [],
                            [RunPar.ExtraMaps[ID].name],
                            [],
                            []
                        ]

    for ID in ramID:
        if ID == 0:
            for reference in FILES.OtherRefs:
                if "RamanAmide" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][2].append("AmideBB")
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            ["AmideBB"],
                            []
                        ]
        else:
            for entry in RunPar.ExtraMaps[ID].references:
                if "Raman" in entry[1]:
                    refstr = str(entry[0])
                    if refstr in allrefs:
                        allrefs[refstr][2].append(RunPar.ExtraMaps[ID].name)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [RunPar.ExtraMaps[ID].name],
                            []
                        ]

    for ID in coupID:
        if ID in [2, 3, 4]:
            coupled = []
            if 0 in freqID:
                coupled.append("AmideBB")
            if 1 in freqID:
                coupled.append(RunPar.ExtraMaps[1].name)
            toapp = []
            for i in coupled:
                for j in coupled:
                    toapp.append(i + " and " + j)
        elif ID in [101, 102]:
            toapp = ["AmideBB and AmideBB"]

        if ID in [0, 1]:
            pass
        elif ID == 2:
            for reference in FILES.OtherRefs:
                if "CoupTDCKrimm" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][3].extend(toapp)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [],
                            toapp
                        ]
        elif ID == 3:
            for reference in FILES.OtherRefs:
                if "CoupTDCTasumi" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][3].extend(toapp)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [],
                            toapp
                        ]
        elif ID == 4:
            for reference in FILES.OtherRefs:
                if "CoupTCC" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][3].extend(toapp)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [],
                            toapp
                        ]
        elif ID == 101:
            for reference in FILES.OtherRefs:
                if "CoupTasumi" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][3].extend(toapp)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [],
                            toapp
                        ]
        elif ID == 102:
            for reference in FILES.OtherRefs:
                if "CoupGLDP" in reference.AIMnotes:
                    refstr = str(reference)
                    if refstr in allrefs:
                        allrefs[refstr][3].extend(toapp)
                    else:
                        allrefs[refstr] = [
                            [],
                            [],
                            [],
                            toapp
                        ]

        else:
            coupled = []
            for pair in RunPar.CouplingMaps[ID].coupled:
                if all(oscID in freqID for oscID in pair):
                    coupled.append(pair)
            for reference in RunPar.CouplingMaps[ID].references:
                refstr = str(reference)
                if refstr in allrefs:
                    allrefs[refstr][3].extend(coupled)
                else:
                    allrefs[refstr] = [
                        [],
                        [],
                        [],
                        coupled
                    ]

    vprint(1, "\nDevelopment of the AIM program:", FILES.logfilename, RunPar)
    for ref in AIMrefs:
        vprint(1, ref, FILES.logfilename, RunPar)


    for refstr, used in allrefs.items():
        vprint(1, "\nFor calculating the:", FILES.logfilename, RunPar)
        if len(used[0]) > 0:
            vprint(
                1, "Frequency of " + ", ".join(used[0]),
                FILES.logfilename, RunPar)
        if len(used[1]) > 0:
            vprint(
                1, "Dipole moment of " + ", ".join(used[1]),
                FILES.logfilename, RunPar
            )
        if len(used[2]) > 0:
            vprint(
                1, "Raman tensor of " + ", ".join(used[2]),
                FILES.logfilename, RunPar
            )
        if len(used[3]) > 0:
            vprint(
                1, "Coupling between " + ", ".join(used[3]),
                FILES.logfilename, RunPar
            )
        vprint(1, refstr, FILES.logfilename, RunPar)
