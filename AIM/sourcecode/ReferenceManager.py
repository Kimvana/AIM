

class Reference():
    def __init__(self, rawref):
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
        # if rawref.startswith("@"):
        #     rawref = rawref[1:]
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
        if len(self.author) > 5:
            return ", ".join(self.author[:5]) + " et al"
        else:
            return ", ".join(self.author)


def readreffile(FILES):
    referencefile = open(FILES.referencefilefilename)
    alldata = referencefile.read()
    allrefs = readrefstring(alldata)
    return allrefs


def readrefstring(string):
    allrefs = []
    temp = string.split("\n@")
    for ref in temp:
        ref = ref.strip()
        if len(ref) > 0:
            refparse = Reference(ref)
            if refparse.success:
                allrefs.append(refparse)

    return allrefs
