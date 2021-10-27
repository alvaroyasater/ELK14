import numpy as np
import pandas as pd
import math as math

# This script computes the load flow with Newton Raphson.
# Buses and lines are organised in classes, and input data is read from two excel files.

#---------------------------------------------------------------------------------------#
# Change these indicators to 1 to include impact of capacitor bank or generator limits.
# In case of an outage on line 2-3, let outage_line_23 = 1 in addition to editing the excel file. This will correct the calculation of line losses.
capacitor_bank = 1
generator_limits = 1
outage_line_23 = 0
#---------------------------------------------------------------------------------------#

# Function to find the complex conjugate
def complex_conj(num):
    conj = np.complex(num.real,-num.imag)
    return conj

class Bus:
    def __init__(self, rowdata):
        self.i = int(rowdata[0]) - 1            # Bus number
        self.pspec = rowdata[1]                 # Specified active load(negative sign) or generation(positive sign)
        self.qspec = rowdata[2]                 # Specified reactive load(negative sign) or generation(positive sign)
        self.v = rowdata[3]                     # Bus voltage magnitude
        self.theta = rowdata[4]                 # Bus voltage angle
        self.type = rowdata[5]                  # Type; slack(type 1), pv(type 2), pq(type 3)
        self.plim = rowdata[6]                  # Active power limit at generator
        self.qlim = rowdata[7]                  # Reactive power limit at generator
        self.shunt = rowdata[8]                 # Shunt capacitor bank added to bus 4 in task A5.

class Line:
    def __init__(self, rowdata):
        self.i = int(rowdata[0]) - 1            # From bus
        self.j = int(rowdata[1]) - 1            # To bus
        self.r = rowdata[2]                     # Line resistance
        self.x = rowdata[3]                     # Line reactance
        self.Yc = 1j * rowdata[4]               # Shunt admittance (also called "line charging capacitance")
        self.y = 1 / (self.r + 1j * self.x)     # Line admittance

#  Read line data from excel file
l = pd.read_excel('linedata.xlsx')
linedata = l.values

#  Create list of Line objects
Lines = []
for row in linedata:
    line = Line(row)
    Lines.append(line)

# Read bus data from excel file
b = pd.read_excel('busdata.xlsx')
busdata = b.values

# Create list of Bus objects
Buses = []
for row in busdata:
    bus = Bus(row)
    Buses.append(bus)

#  Fill Ybus with line (and shunt) admittances. Create matrix of conductance and suseptance.
Ybus = 1j*np.zeros((len(Buses),len(Buses)))
G = np.zeros((len(Buses),len(Buses)))
B = np.zeros((len(Buses),len(Buses)))
np.set_printoptions(precision=2)
for l in Lines:
    Ybus[l.i, l.j] += -l.y
    Ybus[l.j, l.i] += -l.y
    Ybus[l.i, l.i] += l.y + l.Yc / 2
    Ybus[l.j, l.j] += l.y + l.Yc / 2

if capacitor_bank:
    Ybus[3,3] += 1j*Buses[3].shunt      # Add shunt capacitor bank to bus 4 by including its impedance in the Ybus.

for i in range(len(Buses)):
    for j in range(len(Buses)):
        G[i,j] = Ybus[i,j].real
        B[i,j] = Ybus[i,j].imag

# Print for debugging purposes
print('Ybus:')
print(Ybus,'\n')

# Calculating INITIAL power injection for each bus
P = np.zeros((len(Buses)))
Q = np.zeros((len(Buses)))
for i in Buses:
    for j in Buses:
        P[i.i] += i.v * j.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i, j.i] * np.sin(i.theta - j.theta))
        Q[i.i] += i.v * j.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i, j.i] * np.cos(i.theta - j.theta))

# Compute INITIAL change in active and reactive power (dpdq).
# Only include relevant indices (where active and reactive powers are known).
dpdq = []
for i in Buses:
    if i.type == 2 or i.type == 3:
        p = i.pspec - P[i.i]
        dpdq.append(p)
for i in Buses:
    if i.type == 3:
        q = i.qspec - Q[i.i]
        dpdq.append(q)

# Calculating the INITIAL jacobian matrix
# Initialisation of the INITIAL submatirces.
J_PE = np.zeros((len(Buses), len(Buses)))  # dp/dv
J_PT = np.zeros((len(Buses), len(Buses)))  # dp/dtheta
J_QE = np.zeros((len(Buses), len(Buses)))  # dq/dv
J_QT = np.zeros((len(Buses), len(Buses)))  # dq/dtheta

# Construction of the INITIAL submatrices.
for i in Buses:
    # Diagonal elements
    J_PE[i.i, i.i] = P[i.i] / i.v + G[i.i, i.i] * i.v
    J_PT[i.i, i.i] = -Q[i.i] - B[i.i, i.i] * i.v * i.v
    J_QE[i.i, i.i] = Q[i.i] / i.v - B[i.i, i.i] * i.v
    J_QT[i.i, i.i] = P[i.i] - G[i.i, i.i] * i.v * i.v
    for j in Buses:
        if j != i:
            J_PE[i.i, j.i] = i.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i, j.i] * np.sin(i.theta - j.theta))
            J_PT[i.i, j.i] = i.v * j.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i, j.i] * np.cos(i.theta - j.theta))
            J_QE[i.i, j.i] = i.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i, j.i] * np.cos(i.theta - j.theta))
            J_QT[i.i, j.i] = -i.v * j.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i, j.i] * np.sin(i.theta - j.theta))

# Remove extra rows and cols (i.e., unnecessary equations (rows), and variables (cols))
# Rows: slack bus: P & Q, PV buses: Q
# Cols: slack bus: V & theta, PV nodes: V
for i in reversed(Buses):
    if i.type == 2:
        J_QT = np.delete(J_QT, i.i, 0)
        J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)
        J_PE = np.delete(J_PE, i.i, 1)

    if i.type == 1:
        J_PT = np.delete(np.delete(J_PT, i.i, 0), i.i, 1)
        J_PE = np.delete(np.delete(J_PE, i.i, 0), i.i, 1)
        J_QT = np.delete(np.delete(J_QT, i.i, 0), i.i, 1)
        J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)

# Build J as [[J_PT,J_PE],[J_QT,J_QE]]
J_P = np.concatenate((J_PT, J_PE), axis=1)
J_Q = np.concatenate((J_QT, J_QE), axis=1)
J = np.concatenate((J_P, J_Q))


count = 1                                   # Count NR iterations
epsilon = 0                                 # Boolean value to determine convergence. Is later true if convergence criterion is met.
typeswitch = np.zeros((len(Buses)))         # List of boolean values to indicate if a bus is typeswitched (value 1) or not (value 0).
# Start of the iterative NR algorithm.
# Runs until convergence is met or stops at 50 iterations if not convergence before that.
while not epsilon:
    # print('-------------------------------------------------------------------------------')
    # print('Iteration nr %s:' % count)

    # Function to check for Nan in Busdata.xlsx
    def isNaN(num):
        return num != num

    # Compute mismatches in voltage magnitude and angle. x is a list with angles as the first entries and magnitudes at the end.
    x = np.linalg.solve(J,dpdq)
    # print('x: ')
    # print(x)

    # Update new values of theta and v. This will update the objects in the Bus class and in the list Buses.
    k=0
    for i in Buses:
        if i.type == 2 or i.type == 3:
            i.theta += x[k]
            k += 1
    for i in Buses:
        if i.type == 3:
            i.v += x[k]
            k += 1

    #  Calculating power injection for each bus
    P = np.zeros((len(Buses)))
    Q = np.zeros((len(Buses)))
    for i in Buses:
        for j in Buses:
            P[i.i] += i.v * j.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i, j.i] * np.sin(i.theta - j.theta))
            Q[i.i] += i.v * j.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i, j.i] * np.cos(i.theta - j.theta))

    # Check if reactive power limit is exceeded. If yes, typeswitch. If already typeswitched, check if it can be switched back.
    if generator_limits:
        for i in Buses:
            if (typeswitch[i.i] == 0) and (not isNaN(i.qlim)) and Q[i.i] > i.qlim:
                i.qspec = i.qlim
                i.type = 3
                i.v = busdata[i.i, 3]
                typeswitch[i.i] = 1
                print('Bus %s has reached reactive power limit at iteration %d, and is typeswitched to a PQ bus.' % (i.i + 1,count))
            if typeswitch[i.i] == 1 and i.v > busdata[i.i, 3] and Q[i.i] > i.qlim:
                typeswitch[i.i] = 0
                i.qspec = []
                i.v = busdata[i.i, 3]
                print('Bus %s is swiched back to a PV bus.' % (i.i + 1,))

    # Compute change in active and reactive power (dpdq).
    # Only include relevant indices (where active and reactive powers are known).
    dpdq = []
    for i in Buses:
        if i.type == 2 or i.type == 3:
            p = i.pspec -P[i.i]
            dpdq.append(p)
    for i in Buses:
        if i.type == 3:
            q = i.qspec -Q[i.i]
            dpdq.append(q)

    # Calculating the jacobian matrix
    # Initialisation of the submatirces.
    J_PE = np.zeros((len(Buses),len(Buses))) #dp/dv
    J_PT = np.zeros((len(Buses),len(Buses))) #dp/dtheta
    J_QE = np.zeros((len(Buses),len(Buses))) #dq/dv
    J_QT = np.zeros((len(Buses),len(Buses))) #dq/dtheta

    # Construction of the submatrices.
    for i in Buses:
        # Diagonal elements
        J_PE[i.i, i.i] = P[i.i]/i.v + G[i.i,i.i] * i.v
        J_PT[i.i, i.i] = -Q[i.i] - B[i.i, i.i] * i.v * i.v
        J_QE[i.i, i.i] = Q[i.i]/i.v - B[i.i, i.i] * i.v
        J_QT[i.i, i.i] = P[i.i] - G[i.i, i.i] * i.v * i.v
        for j in Buses:
            if j != i:
                J_PE[i.i, j.i] = i.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i,j.i] * np.sin(i.theta - j.theta))
                J_PT[i.i, j.i] = i.v * j.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i, j.i] * np.cos(i.theta - j.theta))
                J_QE[i.i, j.i] = i.v * (G[i.i, j.i] * np.sin(i.theta - j.theta) - B[i.i,j.i] * np.cos(i.theta - j.theta))
                J_QT[i.i, j.i] = -i.v * j.v * (G[i.i, j.i] * np.cos(i.theta - j.theta) + B[i.i, j.i] * np.sin(i.theta - j.theta))

    # Remove extra rows and cols (i.e., unnecessary equations (rows), and variables (cols))
    # Rows: slack bus: P & Q, PV buses: Q
    # Cols: slack bus: V & theta, PV nodes: V
    for i in reversed(Buses):
        if i.type == 2:
            J_QT = np.delete(J_QT, i.i, 0)
            J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)
            J_PE = np.delete(J_PE, i.i, 1)

        if i.type == 1:
            J_PT = np.delete(np.delete(J_PT, i.i, 0), i.i, 1)
            J_PE = np.delete(np.delete(J_PE, i.i, 0), i.i, 1)
            J_QT = np.delete(np.delete(J_QT, i.i, 0), i.i, 1)
            J_QE = np.delete(np.delete(J_QE, i.i, 0), i.i, 1)

    # Build J as [[J_PT,J_PE],[J_QT,J_QE]]
    J_P = np.concatenate((J_PT,J_PE),axis=1)
    J_Q = np.concatenate((J_QT,J_QE),axis=1)
    J = np.concatenate((J_P,J_Q))

    # Check convergence
    epsilon = all(abs(ele) < (10**-5) for ele in dpdq)

    # If convergence, print state of the grid
    if epsilon:
        print('---------------------------------------------------------------------------------------------')
        print('Newton Raphson load flow converges in',count, 'iterations. ')

        # If capacitor bank is connected, include it to the reactive power injection at but 4
        if capacitor_bank:
            Q[3] += -Buses[3].v*Buses[3].v/-1

        if generator_limits:
            print('Generator limits are included.')
        else:
            print('Generator limits are not included.')
        if capacitor_bank:
            print('Capacitor bank is connected at bus 4.')
        else:
            print('Capacitor bank is not connected at bus 4.')
        for i in Buses:
            if typeswitch[i.i] == 1:
                print('Solution reached with bus %s typeswitched to a PQ bus. The bus could not maintain the desired voltage magnitude.\n' % (i.i + 1,))
        print('\nLoad flow solution in pu values:')
        # Print load flow states
        for i in Buses:
            print('| V{:d}'.format(i.i+1), ': {:9f}   '.format(i.v), '| Theta{:d}'.format(i.i+1), ': {:9f}   '.format(i.theta), '| P{:d}'.format(i.i+1), ': {:9f}   '.format(P[i.i]),'| Q{:d}'.format(i.i+1), ': {:9f}  |'.format(Q[i.i]))
        print('\nIn real values:')
        for i in Buses:
            print('| V{:d}'.format(i.i+1), ': {:9f}   '.format(i.v*132), '| Theta{:d}'.format(i.i+1), ': {:9f}   '.format(i.theta), '| P{:d}'.format(i.i+1), ': {:9f}   '.format(P[i.i]*100),'| Q{:d}'.format(i.i+1), ': {:9f}  |'.format(Q[i.i]*100))

        # Compute complex voltages
        v_complex = 1j*np.zeros(len(Buses),dtype=np.complex128)
        v_complex_conj = 1j*np.zeros(len(Buses),dtype=np.complex128)
        for i in Buses:
            v_complex[i.i] = complex(i.v * np.cos(i.theta), i.v * np.sin(i.theta))
            v_complex_conj[i.i] = complex_conj(v_complex[i.i])

        # Compute line flows
        flow = 1j*np.zeros((len(Buses),len(Buses)),dtype=np.complex128)
        flow_oc = 1j * np.zeros((len(Buses), len(Buses)), dtype=np.complex128)
        line_loss_tot = complex(0, 0)
        line_loss_tot_oc = complex(0, 0)
        loss_23_oc = complex(0,0)
        for i in Lines:
            flow[i.i, i.j] = v_complex[i.i] * (v_complex_conj[i.i] - v_complex_conj[i.j]) * (complex_conj(-Ybus[i.i, i.j]))
            flow[i.j, i.i] = v_complex[i.j] * (v_complex_conj[i.j] - v_complex_conj[i.i]) * (complex_conj(-Ybus[i.j, i.i]))
            flow_oc[i.i, i.j] = v_complex[i.i] * (v_complex_conj[i.i] - v_complex_conj[i.j]) * (complex_conj(-Ybus[i.i, i.j])) - v_complex[i.i]*v_complex_conj[i.i]*(i.Yc/2)
            flow_oc[i.j, i.i] = v_complex[i.j] * (v_complex_conj[i.j] - v_complex_conj[i.i]) * (complex_conj(-Ybus[i.j, i.i])) - v_complex[i.j]*v_complex_conj[i.j]*(i.Yc/2)
            line_loss_tot += flow[i.i, i.j] + flow[i.j, i.i]
            line_loss_tot_oc += flow_oc[i.i, i.j] + flow_oc[i.j, i.i]
        if not outage_line_23: #adjusts losses so that line 2-3 is not taken into account twice
            line_loss_tot -= (flow[1, 2] + flow[2, 1])
            line_loss_tot_oc -= (flow_oc[1, 2] + flow_oc[2, 1]) + v_complex[1]*v_complex_conj[1]*(Lines[1].Yc/2) + v_complex[2]*v_complex_conj[2]*(Lines[1].Yc/2)
            flow_oc[1, 2] = v_complex[1] * (v_complex_conj[1] - v_complex_conj[2]) * (complex_conj(-Ybus[1, 2])) - v_complex[1] * v_complex_conj[1] * (Lines[1].Yc)
            flow_oc[2, 1] = v_complex[2] * (v_complex_conj[2] - v_complex_conj[1]) * (complex_conj(-Ybus[2, 1])) - v_complex[2] * v_complex_conj[2] * (Lines[1].Yc)
        loss_23_oc = flow_oc[1, 2] + flow_oc[2, 1]

        # Print complex line power flows
        print('\nLine flow and losses without impact of shunt capacitance:')
        print('| S12:', '{:11f}     '.format(flow[0, 1]), '| S21:', '{:11f}     '.format(flow[1, 0]), '| Loss12:', '{:11f}     '.format(flow[0, 1] + flow[1, 0]))
        print('| S23:', '{:11f}    '.format(flow[1, 2]), '| S32:', '{:11f}      '.format(flow[2, 1]), '| Loss23:', '{:11f}     '.format(flow[1, 2] + flow[2, 1]))
        print('| S14:', '{:11f}     '.format(flow[0, 3]), '| S41:', '{:11f}     '.format(flow[3, 0]), '| Loss14:', '{:11f}     '.format(flow[0, 3] + flow[3, 0]))
        print('| S24:', '{:11f}     '.format(flow[1, 3]), '| S42:', '{:11f}     '.format(flow[3, 1]), '| Loss24:', '{:11f}     '.format(flow[1, 3] + flow[3, 1]))
        print('| S45:', '{:11f}     '.format(flow[3, 4]), '| S54:', '{:11f}     '.format(flow[4, 3]), '| Loss45:', '{:11f}     '.format(flow[3, 4] + flow[4, 3]))
        print('\nLine flow with impact of shunt capacitance:')
        print('| S12:', '{:11f}     '.format(flow_oc[0, 1]), '| S21:', '{:11f}     '.format(flow_oc[1, 0]), '| Loss12:', '{:11f}     '.format(flow_oc[0, 1] + flow_oc[1, 0]))
        print('| S23:', '{:11f}    '.format(flow_oc[1, 2]), '| S32:', '{:11f}      '.format(flow_oc[2, 1]), '| Loss23:', '{:11f}     '.format(loss_23_oc))
        print('| S14:', '{:11f}     '.format(flow_oc[0, 3]), '| S41:', '{:11f}     '.format(flow_oc[3, 0]), '| Loss14:', '{:11f}     '.format(flow_oc[0, 3] + flow_oc[3, 0]))
        print('| S24:', '{:11f}     '.format(flow_oc[1, 3]), '| S42:', '{:11f}     '.format(flow_oc[3, 1]), '| Loss24:', '{:11f}     '.format(flow_oc[1, 3] + flow_oc[3, 1]))
        print('| S45:', '{:11f}     '.format(flow_oc[3, 4]), '| S54:', '{:11f}     '.format(flow_oc[4, 3]), '| Loss45:', '{:11f}     '.format(flow_oc[3, 4] + flow_oc[4, 3]))

        print('\nTotal line losses without operating capacitance: ','{:11f}     '.format(line_loss_tot))
        print('Total line losses with operating capacitance: ', '{:11f}     '.format(line_loss_tot_oc))
        break

    if count > 50:
        print('NR could not converge in 50 iterations...')
        break
    count += 1
    
# SHORT CIRCUIT ANALYSIS


def getVariables(file_shortcircuit):
    #  Read line/bus data from excel file
    df_line_1 = pd.read_excel(file_shortcircuit, 'pos')
    linedata_1 = df_line_1.values
    df_line_2 = pd.read_excel(file_shortcircuit, 'neg')
    linedata_2 = df_line_2.values
    df_line_0 = pd.read_excel(file_shortcircuit, 'zero')
    linedata_0 = df_line_0.values

    #  Create list of TransmissionLine objects
    lineList_1 = []
    for row in linedata_1:
        line = Line(row)
        lineList_1.append(line)

    lineList_2 = []
    for row in linedata_2:
        line = Line(row)
        lineList_2.append(line)

    lineList_0 = []
    for row in linedata_0:
        line = Line(row)
        lineList_0.append(line)

    return lineList_1, lineList_2, lineList_0


def get_zbus012(lineList_1, lineList_2, lineList_0, Buses):

    ybus_1 = 1j * np.zeros((np.size(Buses), np.size(Buses)))
    for k in lineList_1:
        ybus_1[k.i, k.j] -= k.y
        ybus_1[k.j, k.i] -= k.y
        ybus_1[k.i, k.i] += k.y + (k.Yc / 2)
        ybus_1[k.j, k.j] += k.y + (k.Yc / 2)

    ybus_2 = 1j * np.zeros((np.size(Buses), np.size(Buses)))
    for k in lineList_2:
        ybus_2[k.i, k.j] -= k.y
        ybus_2[k.j, k.i] -= k.y
        ybus_2[k.i, k.i] += k.y + (k.Yc / 2)
        ybus_2[k.j, k.j] += k.y + (k.Yc / 2)

    ybus_0 = 1j * np.zeros((np.size(Buses), np.size(Buses)))
    for k in lineList_0:
        ybus_0[k.i, k.j] -= k.y
        ybus_0[k.j, k.i] -= k.y
        ybus_0[k.i, k.i] += k.y + (k.Yc / 2)
        ybus_0[k.j, k.j] += k.y + (k.Yc / 2)
        
    #Adding loads to pos and neg sequence at bus 2, 4 and 5
    ybus_1[1, 1] += (3* (-0.6 + 1j*0.3)) / (Buses[1].v ** 2)
    ybus_1[3, 3] += (3* (-0.6 + 1j*0.2)) / (Buses[3].v ** 2)
    ybus_1[4, 4] += (3* (-0.5 + 1j*0.4)) / (Buses[4].v ** 2)
    ybus_2[1, 1] += (3* (-0.6 + 1j*0.3)) / (Buses[1].v ** 2)
    ybus_2[3, 3] += (3* (-0.6 + 1j*0.2)) / (Buses[3].v ** 2)
    ybus_2[4, 4] += (3* (-0.5 + 1j*0.4)) / (Buses[4].v ** 2)
    
    #Adding generator impedances to all sequences at bus 1 and 3
    ybus_1[0, 0] += 1 / 0.25j
    ybus_2[0, 0] += 1 / 0.15j
    ybus_0[0, 0] += 1 / 0.075j
    ybus_1[2, 2] += 1 / 0.25j
    ybus_2[2, 2] += 1 / 0.15j
    ybus_0[2, 2] += 1 / 0.075j

    #Adding shunt capacitor to all sequences
    ybus_0[3, 3] += 1j
    ybus_1[3, 3] += 1j
    ybus_2[3, 3] += 1j
    
    zbus_1 = np.linalg.inv(ybus_1 )
    zbus_2 = np.linalg.inv(ybus_2)
    zbus_0 = np.linalg.inv(ybus_0)

    return zbus_1, zbus_2, zbus_0


lineList_1, lineList_2, lineList_0 = getVariables('short_circuit.xlsx')
[zbus_1, zbus_2, zbus_0]=get_zbus012(lineList_1, lineList_2, lineList_0, Buses) # as the fault is symmetrical, we only use the positive sequence
print('\n\n SHORT CIRCUIT ANALYSIS: SYMMETRICAL FAULT \n \n Task 7.a)\n\n Positive sequence Zbus:\n')
for i in range(5):
    print('| ', " {:.4f}".format(zbus_1[i,0]), " {:.4f}".format(zbus_1[i,1]), " {:.4f}".format(zbus_1[i,2]), " {:.4f}".format(zbus_1[i,3]), " {:.4f}".format(zbus_1[i,4]), '|')

print('\n\n Task 7.b) \n')
V_prefault_a = []

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

for b in Buses:
    (x,y) = pol2cart(b.v,b.theta)
    V_prefault_a.append(complex(x,y))
    
I_fault_a = V_prefault_a[4]/zbus_1[4, 4]
a = math.cos(2 * math.pi / 3) + 1j * math.sin(2 * math.pi / 3)
I_fault_b = I_fault_a * a**2
I_fault_c = I_fault_a * a
print('Fault currents at bus 5 (in phases a, b and c):\n\n Phase a:  ',  '{:.4} '.format(I_fault_a), '\n Phase b:  ',  '{:.4} '.format(I_fault_b), '\n Phase c:  ', '{:.4} '.format(I_fault_c),'\n\n')

print('Task 7.c) \n')
V_delta_a= np.zeros(5, dtype=complex)
for i in range(5):
    V_delta_a[i] = -(zbus_1[i,4]/zbus_1[4,4])*V_prefault_a[4]

V_postfault_a = V_prefault_a + V_delta_a
V_postfault_b = V_postfault_a * a**2
V_postfault_c = V_postfault_a * a
np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)
x=' '
print('Voltages at all buses (in phases a, b and c):\n \n')
print(18*x, 'Bus1' ,10*x ,'Bus 1' ,10*x ,'Bus 3', 10*x, 'Bus 4', 10*x, 'Bus 5')
print('Phase a:  ', "  {:.4f}".format(V_postfault_a[0]), "  {:.4f}".format(V_postfault_a[1]), "  {:.4f}".format(V_postfault_a[2]), "  {:.4f}".format(V_postfault_a[3]), " {:.4f}".format(V_postfault_a[4]))
print('Phase b:  ', " {:.4f}".format(V_postfault_b[0]), " {:.4f}".format(V_postfault_b[1]), " {:.4f}".format(V_postfault_b[2]), " {:.4f}".format(V_postfault_b[3]), " {:.4f}".format(V_postfault_b[4]))
print('Phase c:  ', " {:.4f}".format(V_postfault_c[0]), " {:.4f}".format(V_postfault_c[1]), " {:.4f}".format(V_postfault_c[2]), " {:.4f}".format(V_postfault_c[3]), " {:.4f}".format(V_postfault_c[4]))

print('\n\nTask 7d)\n')
I_45_a = ((V_postfault_a[3]-V_postfault_a[4])/zbus_1[3,4])
I_45_b = I_45_a * a**2
I_45_c = I_45_a * a
print('Current flowing from bus 4 to 5 in phase a, b, c: \nPhase a: ', " {:.4f}".format(I_45_a), '\nPhase b: ', " {:.4f}".format(I_45_b), '\nPhase c: ', " {:.4f}".format(I_45_c),'\n\n')

print('TASK 8, UNSYMMETRICAL FAULT')

A=np.array([[1,1,1],
   [1,a**2,a],
   [1,a,a**2]], dtype=complex)
A_inv=np.linalg.inv(A)

#8a, negative and zero zequence zbus
print('Task 8a)\n')
[zbus_1, zbus_2, zbus_0]=get_zbus012(lineList_1, lineList_2, lineList_0, Buses)
print('\nPositive sequence Zbus:\n')
for i in range(5):
    print('| ', "{:.4f} ".format(zbus_1[i,0]), "{:.4f} ".format(zbus_1[i,1]), "{:.4f} ".format(zbus_1[i,2]), "{:.4f} ".format(zbus_1[i,3]), "{:.4f} ".format(zbus_1[i,4]), ' |')
print('\nNegative sequence Zbus:\n')
for i in range(5):
    print('| ', "{:.4f} ".format(zbus_2[i,0]), "{:.4f} ".format(zbus_2[i,1]), "{:.4f} ".format(zbus_2[i,2]), "{:.4f} ".format(zbus_2[i,3]), "{:.4f} ".format(zbus_2[i,4]), ' |')
print('\nZero sequence Zbus:\n')
for i in range(5):
    print('| ', "{:.4f} ".format(zbus_0[i,0]), "{:.4f} ".format(zbus_0[i,1]), "{:.4f} ".format(zbus_0[i,2]), "{:.4f} ".format(zbus_0[i,3]), "{:.4f} ".format(zbus_0[i,4]), ' |')
    
#8b, fault currents at bus 5, in phases a,b,c
print('\n\nTask 8b)\n')
print('Fault conditions: \n Ib=0 \n Ic=0 \n Va=Ia*Zf \n')
Z_f_SLG=0
I_b_SLG=0
I_c_SLG=0
I_0_SLG= V_prefault_a[4] / (zbus_0[4,4]+zbus_1[4,4]+zbus_2[4,4]+3*Z_f_SLG) 
I_a_SLG=3*I_0_SLG
I_fault_ph_SLG=[I_a_SLG,0.0,0.0] #Fault currents in phase domain
print('Fault currents in phase a, b, c \n\nPhase a:', " {:.4f}".format(I_a_SLG), '\nPhase b:      ', I_b_SLG, '\nPhase c:      ', I_c_SLG, ' \n\n')

print('Task 8c)\n') #Finding voltages at all buses post fault

print('Postfault voltages at buses: \n')
x=' '
print('Phase:', 7*x, 'a' ,13*x ,'b' ,13*x ,'c')
V_postfault_012=[]
V_postfault_abc=[]
for i in range(5):
    V_postfault_0 = (-zbus_0[i,4]*V_prefault_a[4])/zbus_0[4,4]
    V_postfault_1 = V_prefault_a[4] - (-zbus_1[i,4]*V_prefault_a[4])/zbus_1[4,4]
    V_postfault_2 = (-zbus_2[i,4]*V_prefault_a[4])/zbus_2[4,4]
    V_pf_012 = np.array([V_postfault_0, V_postfault_1, V_postfault_2])
    V_postfault_012.append(V_pf_012)
    V_pf_abc = np.matmul(A,V_pf_012)
    V_postfault_abc.append(V_pf_abc)
    print('Bus', i+1, ':', V_pf_abc)

#Task 8d) current flowing 4 to 5
print('\n\nTask 8d) \n\n')
I_45_0_SLG = (V_postfault_012[3][0]-V_postfault_012[4][0])/zbus_0[3,4]
I_45_1_SLG = (V_postfault_012[3][1]-V_postfault_012[4][1])/zbus_1[3,4]
I_45_2_SLG = (V_postfault_012[3][2]-V_postfault_012[4][2])/zbus_2[3,4]
I_45_012 = np.array([I_45_0_SLG, I_45_1_SLG, I_45_2_SLG])
I_45_abc_SLG = np.matmul(A,I_45_012)

print('Current flowing from bus 4 to 5 is: \n')
print(7*x, 'a' ,12*x ,'b' ,12*x ,'c')
print(I_45_abc_SLG)