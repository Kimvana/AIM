# standard lib imports
import time


class Timer:
    """
    Deals with all timing-related matters. Records and stores when certain
    things happen (which are relative to a weird moment, and therefore not very
    usable). Also calculates the difference between these points to determine
    how long different parts of the calculation took.
    """
    def __init__(self):
        """
        Record when the timer object was created
        """
        self.zero = time.perf_counter()

    def AfterInit(self):
        """
        Record when the initialization phase (the non-repeating part) finished
        """
        self.afterinit = time.perf_counter()
        self.afterinit0 = self.afterinit - self.zero

    def StartCalc(self):
        """
        Record when the difficult part of the code starts
        """
        self.startheavycalc = time.perf_counter() - self.zero

    def Update(self):
        """
        Update the timer (used to see how much time all frames until that point
        have taken)
        """
        self.update = time.perf_counter() - self.zero

    def Endtime(self):
        """
        Record when the calculation finishes.
        """
        self.endtime = time.perf_counter()
        self.endtime0 = self.endtime - self.zero
        self.totheavy = self.endtime0 - self.startheavycalc
