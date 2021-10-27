# ---------------------------------------------------------------------------------------#
# The objects in the system is defined here, as desgnated classes for an object oriented approach.
# More objects can be added, they then have to be added here as well as in the Excel-file.
# ---------------------------------------------------------------------------------------#

class Bus:
    """Class for all bus objects in the system."""
    busCount = 0

    def __init__(self, rowdata, G=0, B=0):
        self.i = int(rowdata[0]) - 1            # Bus number, -1 for indexing purposes
        self.busType = int(rowdata[1])               # Bus Type(1=slack, 2=PV, 3=PQ)
        self.pi_spec = rowdata[2]               # Specified real power injection
        self.qi_spec = rowdata[3]               # Specified imaginary power injection
        self.vmag = rowdata[4]                  # Voltage magnitude
        self.vang = rowdata[5]                  # Voltage angle
        self.G = 0                              # Susceptance Gii, filled in later
        self.B = 0                              # Conductance Bii, filled in later

class Line:
    """Class for all line objects in the system."""
    lineCount = 0

    def __init__(self, rowdata):
        self.i = int(rowdata[0]) - 1            # From bus, -1 for indexing purposes
        self.j = int(rowdata[1]) - 1            # To bus, -1 for indexing purposes
        self.r = rowdata[2]                     # Line resistance
        self.x = rowdata[3]                     # Line reactance
        self.conn = rowdata[4]                  # Indicator (1=connected, 0=disconnected)
        self.y = 1 / complex(self.r, self.x)    # Line admittance
