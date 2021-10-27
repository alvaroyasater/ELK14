#!/usr/bin/python

import numpy as np
from DistribObjects import *
import pandas as pd

from MenuFunctions import ViewFileName

def BuildSystem():
    def renumber(BusList, LineList):
        iloop1 = 0
        sbase = 1   # assume that input values of load are in PU.
        temp = np.zeros(2000,dtype=int)
        while iloop1 < len(BusList):
            obj = BusList[iloop1]
            obj.busext = obj.busnum
            obj.busnum = iloop1 +1
            temp[obj.busext] = obj.busnum
            obj.pload = obj.pload/sbase
            obj.qload = obj.qload / sbase
            iloop1 += 1

        iloop1 = 0
        while iloop1 < len(LineList):
            obj = LineList[iloop1]
            obj.fbus = temp[obj.fbus]
            obj.tbus = temp[obj.tbus]
            iloop1 += 1
        return


    BusList = []
    LineList = []
    file = ViewFileName(filext="xls")
    xls = pd.ExcelFile(file)
    df2 = pd.read_excel(xls, 'Bus')
    values = df2.values
    # Read Bus data  --------------------------------------------
    iloop = 0
    # print(' ')
    while iloop < len(values):
        BusList.append(Bus(busnum=int(values[iloop, 0]), pload=values[iloop, 2], qload=values[iloop, 3]))
        iloop += 1
    df2 = pd.read_excel(xls, 'Branch')
    values = df2.values
    # Read Bus data  --------------------------------------------
    iloop = 0
    # print(' ')
    while iloop < len(values):
        LineList.append(Line(fbus=int(values[iloop, 0]), tbus=int(values[iloop, 1]), r=values[iloop, 2], x=values[iloop, 3], ratea=values[iloop,5], ibstat=int(values[iloop,10])))
        iloop += 1

    renumber(BusList, LineList)

    return BusList, LineList
