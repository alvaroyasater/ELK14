import numpy as np
from System_Objects import *
import pandas as pd
from Menu import ViewFileName

# ---------------------------------------------------------------------------------------#
# This script reads the system data from an Excel-file, and builds the objects for the system.
# ---------------------------------------------------------------------------------------#


def System_Setup():
    """Function to build the system as described in the Excel-file. Returns a list of buses and lines."""

    # Read data from excel file
    file = ViewFileName(filext="xlsx")
    xls = pd.ExcelFile(file)
    bval = pd.read_excel(xls, 'Bus')
    bval = bval.values
    lval = pd.read_excel(xls, 'Branch')
    lval = lval.values

    # Crate lists for buses and branches
    Buses = []
    Lines = []

    # Fill in information from the Excel-file
    for row in bval:
        bus = Bus(row)
        Buses.append(bus)

    for row in lval:
        line = Line(row)
        Lines.append(line)

    return Buses, Lines
