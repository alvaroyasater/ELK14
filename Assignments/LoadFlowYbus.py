import numpy as np
import pandas as pd
import openpyxl
from pandas import DataFrame

from System_Setup import *

# This code performs the Newton-Raphson algorithm on a sample system. The system data is
# read from the Excel-file chosen.

# SYSTEM DECLARATIONS: CHANGE THESE ACCORDING TO THE SYSTEM
# ---------------------------------------------------------------------------------------#
# Any additional lines or buses must be added through the Excel-file.
# A flat start is assumed, otherwise the starting values for voltage magnitudes and angles
# must be specified in the Excel-file.
# The bus-types must also be set in the Excel-file. 1=Slack, 2=PV, 3=PQ.
# Connected lines must be specified by 1 in the "Connect" column in the Excel-file.
# ---------------------------------------------------------------------------------------#

class LoadFlowYbus:
    """Class for calculating the Load Flow solution using the Newton-Raphson method."""

    # Constructor
    def __init__(self, Buses, Lines):
        self.Buses = Buses
        self.Lines = Lines
        self.Pcalc = np.zeros((len(Buses)))
        self.Qcalc = np.zeros((len(Buses)))
        self.vmag = np.ones(len(Buses))
        self.vang = np.zeros(len(Buses))

    # Set up the Y-matrix
    def Ybus(self):
        Ybus = 1j * np.zeros((len(self.Buses), len(self.Buses)))
        for l in self.Lines:
            Ybus[l.i, l.j] += -l.y
            Ybus[l.j, l.i] += -l.y
            Ybus[l.i, l.i] += l.y
            Ybus[l.j, l.j] += l.y
        return Ybus

    # Get the conductance and susceptance matrices
    def GB(self):
        Ybus = self.Ybus()
        G = np.zeros((len(self.Buses), len(self.Buses)))
        B = np.zeros((len(self.Buses), len(self.Buses)))
        for i in range(len(self.Buses)):
            for j in range(len(self.Buses)):
                G[i, j] = Ybus[i, j].real
                B[i, j] = Ybus[i, j].imag
        return G, B

    # Calculating the power injections of each bus
    def PowerInjections(self):
        G, B = self.GB()
        Pcalc = np.zeros((len(self.Buses)))
        Qcalc = np.zeros((len(self.Buses)))
        for i in self.Buses:
            for j in self.Buses:
                Pcalc[i.i] += i.vmag * j.vmag * (
                        G[i.i, j.i] * np.cos(i.vang - j.vang) + B[i.i, j.i] * np.sin(i.vang - j.vang))
                Qcalc[i.i] += i.vmag * j.vmag * (
                        G[i.i, j.i] * np.sin(i.vang - j.vang) - B[i.i, j.i] * np.cos(i.vang - j.vang))
        return Pcalc, Qcalc

    # Print power injections, formatted
    def printPowers(self, powers):
        k = 0
        print('\nCalculated Power Injections:')
        for i in self.Buses:
            if i.busType == 2 or i.busType == 3:
                print('| P{:d}'.format(i.i + 1), ': {:9f}   |'.format(powers[0][k]))
                k += 1
        j = 0
        for i in self.Buses:
            if i.busType == 3:
                print('| Q{:d}'.format(i.i + 1), ': {:9f}   |'.format(powers[1][j]))
                j += 1

    # Get the mismatches in power injections
    def mismatch(self):
        dpdq = []
        for i in self.Buses:
            if i.busType == 2 or i.busType == 3:
                dpdq.append(i.pi_spec - self.Pcalc[i.i])
        for i in self.Buses:
            if i.busType == 3:
                dpdq.append(i.qi_spec - self.Qcalc[i.i])
        return dpdq

    # Print Mismatch Vector, formatted
    def printMismatch(self, dpdq):
        k = 0
        print('\nMismatch Vector for power injections:')
        for i in self.Buses:
            if i.busType == 2 or i.busType == 3:
                print('| dP{:d}'.format(i.i + 1), ': {:9f}   |'.format(dpdq[k]))
                k += 1
        for i in self.Buses:
            if i.busType == 3:
                print('| dQ{:d}'.format(i.i + 1), ': {:9f}   |'.format(dpdq[k]))
                k += 1

    # Calculating the Jacobian matrix
    def Jacobian(self):
        G, B = self.GB()
        # Initialisation of the submatirces.
        J_PE = np.zeros((len(self.Buses), len(self.Buses)))  # dp/dv
        J_PT = np.zeros((len(self.Buses), len(self.Buses)))  # dp/dtheta
        J_QE = np.zeros((len(self.Buses), len(self.Buses)))  # dq/dv
        J_QT = np.zeros((len(self.Buses), len(self.Buses)))  # dq/dtheta

        # Construction of the submatrices.
        for i in self.Buses:
            # Diagonal elements
            J_PE[i.i, i.i] = self.Pcalc[i.i] / i.vmag + G[i.i, i.i] * i.vmag
            J_PT[i.i, i.i] = -self.Qcalc[i.i] - B[i.i, i.i] * i.vmag * i.vmag
            J_QE[i.i, i.i] = self.Qcalc[i.i] / i.vmag - B[i.i, i.i] * i.vmag
            J_QT[i.i, i.i] = self.Pcalc[i.i] - G[i.i, i.i] * i.vmag * i.vmag
            for j in self.Buses:
                if j != i:
                    J_PE[i.i, j.i] = i.vmag * (
                            G[i.i, j.i] * np.cos(i.vang - j.vang) + B[i.i, j.i] * np.sin(i.vang - j.vang))
                    J_PT[i.i, j.i] = i.vmag * j.vmag * (
                            G[i.i, j.i] * np.sin(i.vang - j.vang) - B[i.i, j.i] * np.cos(i.vang - j.vang))
                    J_QE[i.i, j.i] = i.vmag * (
                            G[i.i, j.i] * np.sin(i.vang - j.vang) - B[i.i, j.i] * np.cos(i.vang - j.vang))
                    J_QT[i.i, j.i] = -i.vmag * j.vmag * (
                            G[i.i, j.i] * np.cos(i.vang - j.vang) + B[i.i, j.i] * np.sin(i.vang - j.vang))

        # Remove extra rows and cols (i.e., unnecessary equations (rows), and variables (cols))
        # Rows: slack bus: P & Q, PV buses: Q
        # Cols: slack bus: V & theta, PV nodes: V
        for i in reversed(self.Buses):
            if i.busType == 2:
                J_QT = np.delete(J_QT, i.i, 0)
                J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)
                J_PE = np.delete(J_PE, i.i, 1)

            if i.busType == 1:
                J_PT = np.delete(np.delete(J_PT, i.i, 0), i.i, 1)
                J_PE = np.delete(np.delete(J_PE, i.i, 0), i.i, 1)
                J_QT = np.delete(np.delete(J_QT, i.i, 0), i.i, 1)
                J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)

        # Build J as [[J_PT,J_PE],[J_QT,J_QE]]
        J_P = np.concatenate((J_PT, J_PE), axis=1)
        J_Q = np.concatenate((J_QT, J_QE), axis=1)
        J = np.concatenate((J_P, J_Q))
        return J

    # Help-unit for tij
    def tij(self, gij, bij, ti, tj):
        return gij * np.cos(ti - tj) + bij * np.sin(ti - tj)

    # Help-unit for uij
    def uij(self, gij, bij, ti, tj):
        return gij * np.sin(ti - tj) - bij * np.cos(ti - tj)

    # Print the Jacobian, formatted
    def printJacobian(self, J):
        print('\nJacobian Matrix: ')
        print(DataFrame(J).to_string(index=False, header=False))

    # Get the correction vector for angles and voltages
    def correction(self, J, dpdq):
        x = np.linalg.solve(J, dpdq)
        return x

    # Print the correction vector, formatted
    def printCorrection(self, x):
        print('\nCorrection Vector: ')
        k = 0
        for i in self.Buses:
            if i.busType == 2 or i.busType == 3:
                print('| dT{:d}'.format(i.i + 1), ': {:9f}   |'.format(x[k]))
                k += 1
        for i in self.Buses:
            if i.busType == 3:
                print('| dV{:d}'.format(i.i + 1), ': {:9f}   |'.format(x[k]))
                k += 1

    # Perform one iteration of the NR-algorithm
    def iteration(self):
        J = self.Jacobian()
        self.printJacobian(J)
        powers = self.PowerInjections()
        self.Pcalc = powers[0]
        self.Qcalc = powers[1]
        self.printPowers(powers)
        mismatch = self.mismatch()
        self.printMismatch(mismatch)
        corr = self.correction(J, mismatch)
        self.printCorrection(corr)
        k = 0
        for i in self.Buses:
            if i.busType == 2 or i.busType == 3:
                i.vang += corr[k]
                k += 1
        for i in self.Buses:
            if i.busType == 3:
                i.vmag += corr[k]
                k += 1
        return mismatch

    # Perform the Newton-Raphson iterations until convergence is met
    # Can specify either an iteration limit or a maximum error, or both
    # Default maximum iteration limit is set as 15, default error limit as 10**-5
    def NR(self, itLimit=None, errorLimit=None):
        epsilon = 0
        itCount = 1
        if itLimit is None:
            itLimit = 15
        if errorLimit is None:
            errorLimit = 10**-5

        print('Initial values: ')
        powers = self.PowerInjections()
        self.Pcalc = powers[0]
        self.Qcalc = powers[1]
        self.printPowers(powers)
        mismatch = self.mismatch()
        self.printMismatch(mismatch)

        while not epsilon:

            print('-------------------------------------------------------------------------------------------')
            print('Iteration Number: ', itCount)
            mismatch = self.iteration()

            # Check convergence
            epsilon = all(abs(ele) < errorLimit for ele in mismatch)

            if epsilon or itCount >= itLimit:
                if epsilon:
                    print('-------------------------------------------------------------------------------------------')
                    print('Newton Raphson load flow converges in', itCount, 'iterations. ')
                if itCount >= itLimit:
                    print('-------------------------------------------------------------------------------------------')
                    print('Iteration limit of ', itCount, ' iterations was reached')
                print('\nLoad flow solution in pu values:')
                # Print load flow states
                for i in self.Buses:
                    print('| V{:d}'.format(i.i + 1), ': {:9f}   '.format(i.vmag), '| Theta{:d}'.format(i.i + 1),
                          ': {:9f}   '.format(i.vang), '| P{:d}'.format(i.i + 1), ': {:9f}   '.format(self.Pcalc[i.i]),
                          '| Q{:d}'.format(i.i + 1), ': {:9f}  |'.format(self.Qcalc[i.i]))
                break
            itCount += 1


Buses, Lines = System_Setup()
LF = LoadFlowYbus(Buses, Lines)
LF.NR(itLimit=4, errorLimit=10**-5)
