# Copyright (c) 2021, Olav B. Fosso, NTNU
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright notice,
#       this list of conditions and the following disclaimer in the documentation
#       and/or other materials provided with the distribution.
import math

import matplotlib.pyplot as plt

import ElkObjects
from ElkBuildSystem import *


class LoadFlow:
    """
    Common base class   Load Flow
    Input:
        BusList      - List off all Bus objects
        LineList    - List of all transmission lines objects
    Returns: None
    """

    def __init__(self, Buses, Lines):
        self.BusList = Buses
        self.LineList = Lines
        self.voang = np.zeros(len(self.BusList))
        self.vomag = np.ones(len(self.BusList))
        self.topology = []


# Add methods as needed






# -------------------------------------------------------------------------------------------------------------
#
# Demo case (Illustration of how to build up a script)
#
BusList, LineList = BuildSystem()  # Import data from Excel file
# BusList and LineList contain lists of Bus-objects and Line-objects
lf = LoadFlow(BusList, LineList)  # Create object

# How to access objects
for x in BusList:
    print('Busname :', x.busname)

for x in LineList:
    print('Frombus :', x.fbus, '  ToBus :', x.tbus)
