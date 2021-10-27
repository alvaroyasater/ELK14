import numpy as np
from pandas import DataFrame

from System_Setup import System_Setup

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

class LoadFlowLineDriven:
    """Class for calculating the Load Flow solution using the Newton-Raphson method."""

    # Constructor
    def __init__(self, Buses, Lines):
        self.Buses = Buses
        self.Lines = Lines
        self.Pcalc = np.zeros((len(Buses)))
        self.Qcalc = np.zeros((len(Buses)))
        self.vmag = np.ones(len(Buses))
        self.vang = np.zeros(len(Buses))
        self.GiiBii()
        self.JacobianLine()

    # Finding Gii and Bii for all buses
    def GiiBii(self):
        for l in self.Lines:
            self.Buses[l.i].G += l.y.real
            self.Buses[l.j].G += l.y.real
            self.Buses[l.i].B += l.y.imag
            self.Buses[l.j].B += l.y.imag

    # Calculating power injections without Y-bus
    def PowerInjectionsLine(self):
        Pcalc = np.zeros((len(self.Buses)))
        Qcalc = np.zeros((len(self.Buses)))
        for b in self.Buses:
            Pcalc[b.i] = b.vmag * b.vmag * b.G
            Qcalc[b.i] = -b.vmag * b.vmag * b.B
            for l in self.Lines:
                if b.i == l.j:
                    Pcalc[b.i] -= b.vmag * self.Buses[l.i].vmag * self.tij(l.y.real, l.y.imag, b.vang, self.Buses[l.i].vang)
                    Qcalc[b.i] -= b.vmag * self.Buses[l.i].vmag * self.uij(l.y.real, l.y.imag, b.vang, self.Buses[l.i].vang)
                elif b.i == l.i:
                    Pcalc[b.i] -= b.vmag * self.Buses[l.j].vmag * self.tij(l.y.real, l.y.imag, b.vang, self.Buses[l.j].vang)
                    Qcalc[b.i] -= b.vmag * self.Buses[l.j].vmag * self.uij(l.y.real, l.y.imag, b.vang, self.Buses[l.j].vang)
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

    # Help-unit for tij
    def tij(self, gij, bij, ti, tj):
        return gij * np.cos(ti - tj) + bij * np.sin(ti - tj)

    # Help-unit for uij
    def uij(self, gij, bij, ti, tj):
        return gij * np.sin(ti - tj) - bij * np.cos(ti - tj)

    # Line driven Jacobian
    def JacobianLine(self):
        # Initialising the submatrices of Jacobian
        J_PT = np.zeros((len(Buses), len(Buses)))  # dp/dtheta
        J_PE = np.zeros((len(Buses), len(Buses)))  # dp/dv
        J_QT = np.zeros((len(Buses), len(Buses)))  # dq/dtheta
        J_QE = np.zeros((len(Buses), len(Buses)))  # dq/dv

        # Go through all lines and fill submatrices of Jacobian
        for l in self.Lines:
            busi = self.Buses[l.i]
            busj = self.Buses[l.j]

            # Define the help units
            uij = self.uij(l.y.real, l.y.imag, busi.vang, busj.vang)
            tij = self.tij(l.y.real, l.y.imag, busi.vang, busj.vang)
            uji = self.uij(l.y.real, l.y.imag, busj.vang, busi.vang)
            tji = self.tij(l.y.real, l.y.imag, busj.vang, busi.vang)

            # Filling in partial derivatives:
            J_PT[l.i][l.i] += busi.vmag * busj.vmag * uij
            J_PT[l.j][l.j] += busi.vmag * busj.vmag * uji
            J_PT[l.i][l.j] = - busi.vmag * busj.vmag * uij
            J_PT[l.j][l.i] = - busi.vmag * busj.vmag * uji

            J_PE[l.i][l.i] -= busj.vmag * tij
            J_PE[l.i][l.j] = - busi.vmag * tij
            J_PE[l.j][l.j] -= busi.vmag * tji
            J_PE[l.j][l.i] = - busj.vmag * tji

            J_QT[l.i][l.i] -= busi.vmag * busj.vmag * tij
            J_QT[l.j][l.j] -= busi.vmag * busj.vmag * tji
            J_QT[l.i][l.j] = busi.vmag * busj.vmag * tij
            J_QT[l.j][l.i] = busi.vmag * busj.vmag * tji

            J_QE[l.i][l.i] -= busj.vmag * uij
            J_QE[l.j][l.j] -= busi.vmag * uji
            J_QE[l.i][l.j] = - busi.vmag * uij
            J_QE[l.j][l.i] = - busj.vmag * uji

        # Add the elements for Gii and Bii
        for i in self.Buses:
            J_PE[i.i][i.i] += 2 * self.Buses[i.i].vmag * self.Buses[i.i].G
            J_QE[i.i][i.i] -= 2 * self.Buses[i.i].vmag * self.Buses[i.i].B

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
        self.Jacobian = J
        return J

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
        J = self.JacobianLine()
        self.printJacobian(J)
        powers = self.PowerInjectionsLine()
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
        powers = self.PowerInjectionsLine()
        self.Pcalc = powers[0]
        self.Qcalc = powers[1]
        self.printPowers(powers)
        mismatch = self.mismatch()
        self.printMismatch(mismatch)

        while not epsilon:
            print('-------------------------------------------------------------------------------------------')
            print('Iteration Number: ', itCount)
            mismatch = self.iteration()

            # Check for convergence
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
LF = LoadFlowLineDriven(Buses, Lines)
LF.NR(itLimit=4, errorLimit=10**-5)
