

class Reference():
    """
    This class stores all information on a reference. Calling the str()
    function on an instance of this class returns a nice printable format.
    While a new instance could be created by calling this class directly with
    a string in .bib format, the class is intended to be used together with the
    readreffile and readrefstring functions in this file.
    """
    def __init__(self, rawref):
        """
        Manages all functions in this class
        """
        self.dirattr = self.getdirattr()
        self.allattr = self.getdirattr() + ["authorraw", "AIMnotes"]

        try:
            self.strtoattr(rawref)
        except Exception:
            pass

        if any(hasattr(self, x) for x in self.allattr):
            self.success = True
        else:
            self.success = False

        if hasattr(self, "authorraw"):
            self.formatauthor()

    def strtoattr(self, rawref):
        """
        Takes the input string and separates it into all separate entries.
        If an entry is required (some are ignored), it is saved as an attribute
        of this class.
        """
        rawreflist = rawref.split("\n")[1:]
        for line in rawreflist:
            linelist = line.split("=")
            linelist = [item.strip("\t {{}},") for item in linelist]
            if linelist[0] in self.dirattr:
                setattr(self, linelist[0], linelist[1])
            elif linelist[0] == "author":
                authorlistraw = linelist[1].split(" and ")
                authorlist = []
                for author in authorlistraw:
                    names = author.split(", ")
                    authorlist.append([names[1], names[0]])
                self.authorraw = authorlist
            elif linelist[0] == "AIMnotes":
                self.AIMnotes = linelist[1].split(" ")

    def getdirattr(self):
        """
        Defines what kind of information should be collected directly from the
        reference string.
        """
        dirattr = [
            "title",
            "year",
            "month",
            "journal",
            "volume",
            "number",
            "pages",
            "issn",
            "doi",
            "url"
        ]
        return dirattr

    def formatauthor(self):
        """
        Takes the list of authors and converts it into a more usable format.
        """
        self.author = []
        for author in self.authorraw:
            first_name = author[0]
            last_name = author[1]
            firstlist = first_name.split(" ")
            templist = []
            for name in firstlist:
                templist.append(name[0] + ".")
            self.author.append(" ".join(templist) + " " + last_name)

    def __str__(self):
        """
        Makes a nicely formatted string for printing. Contains all information
        from this class. For some attributes, if they're missing, that will be
        indicated. For others, missing information will just not be printed.
        """
        outstr = ""
        if hasattr(self, "author"):
            outstr += self.authorstr() + ". "
        else:
            outstr += "unknown author. "

        if hasattr(self, "title"):
            outstr += '"' + self.title + '". '
        else:
            outstr += "unknown title. "

        outstr += "In "
        outstr += getattr(self, "journal", "unknown journal") + " "

        if hasattr(self, "volume"):
            if hasattr(self, "number"):
                outstr += self.volume + "." + self.number
            else:
                outstr += "vol. " + self.volume
        else:
            if hasattr(self, "number"):
                outstr += "no. " + self.number
            else:
                outstr += "unknown volume"

        if hasattr(self, "year"):
            if hasattr(self, "month"):
                outstr += " (" + self.month + " " + self.year + ")"
            else:
                outstr += " (" + self.year + ")"
        else:
            outstr += " unknown year"

        if hasattr(self, "pages"):
            if "-" in self.pages:
                outstr += ", pp. " + self.pages
            else:
                outstr += ", p. " + self.pages

        for type in ["issn", "doi", "url"]:
            if hasattr(self, type):
                outstr += ", " + type + ": " + getattr(self, type)

        return outstr

    def authorstr(self):
        """
        Converts the author list into something usable for printing.
        """
        if len(self.author) > 5:
            return ", ".join(self.author[:5]) + " et al"
        else:
            return ", ".join(self.author)


def readreffile(reffilefilename):
    """
    Takes a filename (/location), reads the contents, and extracts the
    references from it.
    Returns a list of instances of the Reference class.
    """
    referencefile = open(reffilefilename)
    alldata = referencefile.read()
    allrefs = readrefstring(alldata)
    return allrefs


def readrefstring(string):
    """
    Takes a string (may contain multiple references), and extracts the
    references from it. Returns a list of instances of the Reference class.
    """
    allrefs = []
    temp = string.split("\n@")
    for ref in temp:
        ref = ref.strip()
        if len(ref) > 0:
            refparse = Reference(ref)
            if refparse.success:
                allrefs.append(refparse)

    return allrefs
