import numpy as np
import matplotlib.pyplot as plt
import sys


def z(r,x):
    return complex(r,x)
def bij(r,x):
    return (1.0/complex(r,x)).imag
def gij(r,x):
    return (1.0/complex(r,x)).real

def uij(gij,bij,tetai,tetaj):
    return (gij*np.sin(tetai-tetaj)-bij*np.cos(tetai-tetaj))

def tij(gij,bij,tetai,tetaj):
    return (gij*np.cos(tetai-tetaj)+bij*np.sin(tetai-tetaj))

def uji(gij,bij,tetai,tetaj):
    return (gij*np.sin(tetaj-tetai)-bij*np.cos(tetaj-tetai))

def tji(gij,bij,tetai,tetaj):
    return (gij*np.cos(tetaj-tetai)+bij*np.sin(tetaj-tetai))

r12 = 0.1
x12 = 0.2
r13 = 0.05
x13 = 0.25
r23 = 0.05
x23 = 0.15
v1 = 1.0
v2 = 1.0
v3 = 1.0
teta1 = 0.0
teta2 = 0.0
teta3 = 0.0

g12 = gij(r12, x12)
b12 = bij(r12, x12)
g13 = gij(r13, x13)
b13 = bij(r13, x13)
g23 = gij(r23, x23)
b23 = bij(r23, x23)

# Loads
P01 = 0.8
Q01 = 0.5
P02 = 0.4
Q02 = 0.5

# Define useful functions

def PQ_injections(v1,v2,v3,teta1,teta2,teta3):
    p1 = v1 * v1 * (g12 + g13) - v1 * v2 * tij(g12, b12, teta1, teta2) - v1 * v3 * tij(g13, b13, teta1, teta3)
    p2 = v2 * v2 * (g12 + g23) - v1 * v2 * tij(g12, b12, teta2, teta1) - v2 * v3 * tij(g23, b23, teta2, teta3)
    p3 = v3 * v3 * (g23 + g13) - v1 * v3 * tij(g13, b13, teta3, teta1) - v2 * v3 * tij(g23, b23, teta3, teta2)
    q1 = -v1 * v1 * (b12 + b13) - v1 * v2 * uij(g12, b12, teta1, teta2) - v1 * v3 * uij(g13, b13, teta1, teta3)
    q2 = -v2 * v2 * (b12 + b23) - v1 * v2 * uij(g12, b12, teta2, teta1) - v2 * v3 * uij(g23, b23, teta2, teta3)
    q3 = -v3 * v3 * (b23 + b13) - v1 * v3 * uij(g13, b13, teta3, teta1) - v2 * v3 * uij(g23, b23, teta3, teta2)
    return p1, p2, p3, q1, q2, q3

# dP/dTeta
def dPdTeta(v1,v2,v3,teta1,teta2,teta3):
    dp1dt1 = v1 * v2 * uij(g12, b12, teta1, teta2) + v1 * v3 * uij(g13, b13, teta1, teta3)
    dp1dt2 = -v1 * v2 * uij(g12, b12, teta1, teta2)
    dp1dt3 = -v1 * v3 * uij(g13, b13, teta1, teta3)

    dp2dt2 = v1 * v2 * uij(g12, b12, teta2, teta1) + v2 * v3 * uij(g23, b23, teta2, teta3)
    dp2dt1 = -v1 * v2 * uij(g12, b12, teta2, teta1)
    dp2dt3 = -v2 * v3 * uij(g23, b23, teta2, teta3)

    dp3dt3 = v1 * v3 * uij(g13, b13, teta3, teta1) + v2 * v3 * uij(g23, b23, teta3, teta2)
    dp3dt1 = -v1 * v3 * uij(g13, b13, teta3, teta1)
    dp3dt2 = -v2 * v3 * uij(g23, b23, teta3, teta2)
    return dp1dt1, dp1dt2, dp1dt3, dp2dt1, dp2dt2, dp2dt3, dp3dt1, dp3dt2, dp3dt3

# dP/dV
def dPdV(v1,v2,v3,teta1,teta2,teta3):
    dp1dv1 = 2 * v1 * (g12 + g13) - v2 * tij(g12, b12, teta1, teta2) - v3 * tij(g13, b13, teta1, teta3)
    dp1dv2 = -v1 * tij(g12, b12, teta1, teta2)
    dp1dv3 = -v1 * tij(g13, b13, teta1, teta3)

    dp2dv2 = 2 * v2 * (g12 + g23) - v1 * tij(g12, b12, teta2, teta1) - v3 * tij(g23, b23, teta2, teta3)
    dp2dv1 = -v2 * tij(g12, b12, teta2, teta1)
    dp2dv3 = -v2 * tij(g23, b23, teta2, teta3)

    dp3dv3 = 2 * v3 * (g13 + g23) - v1 * tij(g13, b13, teta3, teta1) - v2 * tij(g23, b23, teta3, teta2)
    dp3dv1 = -v3 * tij(g13, b13, teta3, teta1)
    dp3dv2 = -v3 * tij(g23, b23, teta3, teta2)
    return dp1dv1, dp1dv2, dp1dv3, dp2dv1, dp2dv2, dp2dv3, dp3dv1, dp3dv2, dp3dv3

def dQdTeta(v1, v2, v3, teta1, teta2, teta3):
    dq1dt1 = -v1 * v2 * tij(g12, b12, teta1, teta2) - v1 * v3 * tij(g13, b13, teta1, teta3)
    dq1dt2 = v1 * v2 * tij(g12, b12, teta1, teta2)
    dq1dt3 = v1 * v3 * tij(g13, b13, teta1, teta3)

    dq2dt2 = -v1 * v2 * tij(g12, b12, teta2, teta1) - v2 * v3 * tij(g23, b23, teta2, teta3)
    dq2dt1 = v1 * v2 * tij(g12, b12, teta2, teta1)
    dq2dt3 = v1 * v3 * tij(g23, b23, teta2, teta3)

    dq3dt3 = -v1 * v3 * tij(g13, b13, teta3, teta1) - v2 * v3 * tij(g23, b23, teta3, teta2)
    dq3dt1 = v1 * v2 * tij(g13, b13, teta3, teta1)
    dq3dt2 = v2 * v3 * tij(g23, b23, teta3, teta2)
    return dq1dt1, dq1dt2, dq1dt3, dq2dt1, dq2dt2, dq2dt3, dq3dt1, dq3dt2, dq3dt3

def dQdV(v1, v2, v3, teta1, teta2, teta3):
    dq1dv1 = -2 * v1 * (b12 + b13) - v2 * uij(g12, b12, teta1, teta2) - v3 * uij(g13, b13, teta1, teta3)
    dq1dv2 = -v1 * uij(g12, b12, teta1, teta2)
    dq1dv3 = -v1 * uij(g13, b13, teta1, teta3)

    dq2dv2 = -2 * v2 * (b12 + b23) - v1 * uij(g12, b12, teta2, teta1) - v3 * uij(g23, b23, teta2, teta3)
    dq2dv1 = -v2 * uij(g12, b12, teta2, teta1)
    dq2dv3 = -v2 * uij(g23, b23, teta2, teta3)

    dq3dv3 = -2 * v3 * (b13 + b23) - v1 * uij(g13, b13, teta3, teta1) - v2 * uij(g23, b23, teta3, teta2)
    dq3dv1 = -v3 * uij(g13, b13, teta3, teta1)
    dq3dv2 = -v3 * uij(g23, b23, teta3, teta2)
    return dq1dv1, dq1dv2, dq1dv3, dq2dv1, dq2dv2, dq2dv3, dq3dv1, dq3dv2, dq3dv3


# Print partial deriviative (Jacobian)
def print_jacobian_to_screen():
    print('  dp1/dt1 :', '{:9.4f}'.format(dp1dt1), '  dp1/dt2 :','{:9.4f}'.format(dp1dt2), '  dp1/dt3 :', '{:9.4f}'.format(dp1dt3))
    print('  dp2/dt1 :', '{:9.4f}'.format(dp2dt1), '  dp2/dt2 :','{:9.4f}'.format(dp2dt2), '  dp2/dt3 :', '{:9.4f}'.format(dp2dt3))
    print('  dp3/dt1 :', '{:9.4f}'.format(dp3dt1), '  dp3/dt2 :','{:9.4f}'.format(dp3dt2), '  dp3/dt3 :', '{:9.4f}'.format(dp3dt3))
    print(' ')
    print('  dp1/dv1 :', '{:9.4f}'.format(dp1dv1), '  dp1/dv2 :','{:9.4f}'.format(dp1dv2), '  dp1/dv3 :', '{:9.4f}'.format(dp1dv3))
    print('  dp2/dv1 :', '{:9.4f}'.format(dp2dv1), '  dp2/dv2 :','{:9.4f}'.format(dp2dv2), '  dp2/dv3 :', '{:9.4f}'.format(dp2dv3))
    print('  dp3/dv1 :', '{:9.4f}'.format(dp3dv1), '  dp3/dv2 :','{:9.4f}'.format(dp3dv2), '  dp3/dv3 :', '{:9.4f}'.format(dp3dv3))
    print(' ')
    print('  dq1/dt1 :', '{:9.4f}'.format(dq1dt1), '  dq1/dt2 :','{:9.4f}'.format(dq1dt2), '  dq1/dt3 :', '{:9.4f}'.format(dq1dt3))
    print('  dq2/dt1 :', '{:9.4f}'.format(dq2dt1), '  dq2/dt2 :','{:9.4f}'.format(dq2dt2), '  dq2/dt3 :', '{:9.4f}'.format(dq2dt3))
    print('  dq3/dt1 :', '{:9.4f}'.format(dq3dt1), '  dq3/dt2 :','{:9.4f}'.format(dq3dt2), '  dq3/dt3 :', '{:9.4f}'.format(dq3dt3))
    print(' ')
    print('  dq1/dv1 :', '{:9.4f}'.format(dq1dv1), '  dq1/dv2 :','{:9.4f}'.format(dq1dv2), '  dq1/dv3 :', '{:9.4f}'.format(dq1dv3))
    print('  dq2/dv1 :', '{:9.4f}'.format(dq2dv1), '  dq2/dv2 :','{:9.4f}'.format(dq2dv2), '  dq2/dv3 :', '{:9.4f}'.format(dq2dv3))
    print('  dq3/dv1 :', '{:9.4f}'.format(dq3dv1), '  dq3/dv2 :','{:9.4f}'.format(dq3dv2), '  dq3/dv3 :', '{:9.4f}'.format(dq3dv3))
    return

# Write the partial derivatives to file (Jacobian)
def write_jacobian_to_file():
    f.write('\n' + 'Iteration number  : ' + str(i) + '\n\n')
    f.write('  dp1/dt1 :' + '{:9.4f}'.format(dp1dt1) + '  dp1/dt2 :' + '{:9.4f}'.format(
        dp1dt2) + '  dp1/dt3 :' + '{:9.4f}'.format(dp1dt3) + '\n')
    f.write('  dp2/dt1 :' + '{:9.4f}'.format(dp2dt1) + '  dp2/dt2 :' + '{:9.4f}'.format(
        dp2dt2) + '  dp2/dt3 :' + '{:9.4f}'.format(dp2dt3) + '\n')
    f.write('  dp3/dt1 :' + '{:9.4f}'.format(dp3dt1) + '  dp3/dt2 :' + '{:9.4f}'.format(
        dp3dt2) + '  dp3/dt3 :' + '{:9.4f}'.format(dp3dt3) + '\n')
    f.write('\n' + ' Next block ' + '\n\n')
    f.write('  dp1/dv1 :' + '{:9.4f}'.format(dp1dv1) + '  dp1/dv2 :' + '{:9.4f}'.format(
        dp1dv2) + '  dp1/dv3 :' + '{:9.4f}'.format(dp1dv3) + '\n')
    f.write('  dp2/dv1 :' + '{:9.4f}'.format(dp2dv1) + '  dp2/dv2 :' + '{:9.4f}'.format(
        dp2dv2) + '  dp2/dv3 :' + '{:9.4f}'.format(dp2dv3) + '\n')
    f.write('  dp3/dv1 :' + '{:9.4f}'.format(dp3dv1) + '  dp3/dv2 :' + '{:9.4f}'.format(
        dp3dv2) + '  dp3/dv3 :' + '{:9.4f}'.format(dp3dv3) + '\n')
    f.write('\n' + ' Next block' + '\n\n')
    f.write('  dq1/dt1 :' + '{:9.4f}'.format(dq1dt1) + '  dq1/dt2 :' + '{:9.4f}'.format(
        dq1dt2) + '  dq1/dt3 :' + '{:9.4f}'.format(dq1dt3) + '\n')
    f.write('  dq2/dt1 :' + '{:9.4f}'.format(dq2dt1) + '  dq2/dt2 :' + '{:9.4f}'.format(
        dq2dt2) + '  dq2/dt3 :' + '{:9.4f}'.format(dq2dt3) + '\n')
    f.write('  dq3/dt1 :' + '{:9.4f}'.format(dq3dt1) + '  dq3/dt2 :' + '{:9.4f}'.format(
        dq3dt2) + '  dq3/dt3 :' + '{:9.4f}'.format(dq3dt3) + '\n')
    f.write('\n' + ' Next block' + '\n\n')
    f.write('  dq1/dv1 :' + '{:9.4f}'.format(dq1dv1) + '  dq1/dv2 :' + '{:9.4f}'.format(
        dq1dv2) + '  dq1/dv3 :' + '{:9.4f}'.format(dq1dv3) + '\n')
    f.write('  dq2/dv1 :' + '{:9.4f}'.format(dq2dv1) + '  dq2/dv2 :' + '{:9.4f}'.format(
        dq2dv2) + '  dq2/dv3 :' + '{:9.4f}'.format(dq2dv3) + '\n')
    f.write('  dq3/dv1 :' + '{:9.4f}'.format(dq3dv1) + '  dq3/dv2 :' + '{:9.4f}'.format(
        dq3dv2) + '  dq3/dv3 :' + '{:9.4f}'.format(dq3dv3) + '\n')
    return

#  Save to file (in case you system reports error - check the 'opsys'
opsys = sys.platform
if opsys == 'linux':
    f = open("/home/olav/Dropbox/Python_code/elk14/cpf_results.txt","w+")
elif opsys == 'windows':
    f = open("cpf_results_new.txt","w+")
elif opsys == 'win32':
    f = open("cpf_results_new.txt","w+")

i = 0
#----------------------------------------------------------------- Iterative Load Flow Part -------------------------------------
#  Iterative loop for load flow calculations
while i < 4:
    # Calculate net injections
    print('-------------------Iteration ', i, '-------------------')
    p1, p2, p3, q1, q2, q3 = PQ_injections(v1, v2, v3, teta1, teta2, teta3)
    print(' ')
    #   print('Injction: ', pA1,p2,p3,q1,q2,q3)
    print('Net injection   P1 :', '{:9.5f}'.format(p1), '  P2 :','{:9.5f}'.format(p2), '  P3 :', '{:9.5f}'.format(p3), '\n',
          '               Q1 :', '{:9.5f}'.format(q1),  '  Q2 :', '{:9.5f}'.format(q2),  '  Q3 :', '{:9.5f}'.format(q3))
    f.write('\n' + 'Net injection   P1 :'+ '{:9.5f}'.format(p1) + '  P2 :'+'{:9.5f}'.format(p2) + '  P3 :'+ '{:9.5f}'.format(p3)+ '\n'+
            '                Q1 :'+ '{:9.5f}'.format(q1) +  '  Q2 :'+ '{:9.5f}'.format(q2)+  '  Q3 :'+ '{:9.5f}'.format(q3) +'\n\n')
    print(' ')

    #Jacobi elements
    #Active power
    dp1dt1, dp1dt2, dp1dt3, dp2dt1, dp2dt2, dp2dt3, dp3dt1, dp3dt2, dp3dt3 = dPdTeta(v1, v2, v3, teta1, teta2, teta3)
    dp1dv1, dp1dv2, dp1dv3, dp2dv1, dp2dv2, dp2dv3, dp3dv1, dp3dv2, dp3dv3 = dPdV(v1, v2, v3, teta1, teta2, teta3)

    # Reactive power
    dq1dt1, dq1dt2, dq1dt3, dq2dt1, dq2dt2, dq2dt3, dq3dt1, dq3dt2, dq3dt3 = dQdTeta(v1, v2, v3, teta1, teta2, teta3)
    dq1dv1, dq1dv2, dq1dv3, dq2dv1, dq2dv2, dq2dv3, dq3dv1, dq3dv2, dq3dv3 = dQdV(v1, v2, v3, teta1, teta2, teta3)



    # Display the Jacobian elements (screen and file
    print_jacobian_to_screen()
    write_jacobian_to_file()

    jac = np.array(
        [[dp1dt1, dp1dt2, dp1dv1, dp1dv2], [dp2dt1, dp2dt2, dp2dv1, dp2dv2], [dq1dt1, dq1dt2, dq1dv1, dq1dv2],
         [dq2dt1, dq2dt2, dq2dv1, dq2dv2]])
    b = np.array([-p1-P01, -p2-P02, -q1-Q01, -q2-Q02])
    #    print(' ')
    #    print('b-vector: ', b)
    print('RHS-vector   dp1 :', '{:10.7f}'.format(-p1-P01), '  dp2 :','{:10.7f}'.format(-p2-P02), '\n',
          '            dq1 :', '{:10.7f}'.format(-q1-Q01),  '  dq2 :', '{:10.7f}'.format(-q2-Q02))
    f.write('\n' +'RHS-vector   dp1 :'+ '{:10.7f}'.format(-p1-P01)+ '  dp2 :'+'{:10.7f}'.format(-p2-P02)+ '\n'+
            '             dq1 :'+ '{:10.7f}'.format(-q1-Q01)+  '  dq2 :'+ '{:10.7f}'.format(-q2-Q02) + '\n\n')
    print(' ')

    x = np.linalg.solve(jac,b)

    #   print(' ')
    #    print('x-vector: ', x)
    print('Correction   dt1 :', '{:10.7f}'.format(x[0]), '  dt2 :','{:10.7f}'.format(x[1]), '\n',
          
          
          '            dv1 :', '{:10.7f}'.format(x[2]),  '  dv2 :', '{:10.7f}'.format(x[3]))
    f.write('\n' + 'Correction   dt1 :'+ '{:10.7f}'.format(x[0])+ '  dt2 :'+'{:10.7f}'.format(x[1])+ '\n'+
            '             dv1 :'+ '{:10.7f}'.format(x[2])+  '  dv2 :'+ '{:10.7f}'.format(x[3]) + '\n\n')
    print(' ')

    # Update the voltage angles and magnitudes
    teta1 = teta1 + x[0]
    teta2 = teta2 + x[1]
    v1 = v1 + x[2]
    v2 = v2 + x[3]
    i = i + 1

# ------------------------- Start the Continuation Power Flow part -------------------------------------
# The Jacobian may be used from the previous iteration

contpar = 'volt'

# Define load distribution
beta = np.array([0.5, 0.5])
alfa = np.array([0.0, 0.0])

jac = np.array(
    [[dp1dt1, dp1dt2, dp1dv1, dp1dv2, beta[0]], [dp2dt1, dp2dt2, dp2dv1, dp2dv2, beta[1]], [dq1dt1, dq1dt2, dq1dv1, dq1dv2, alfa[0]],
     [dq2dt1, dq2dt2, dq2dv1, dq2dv2, alfa[1]], [0.0, 0.0, 0.0, 0.0, 1.0]])
b = np.array([0.0, 0.0, 0.0, 0.0, 1.0])

x = np.linalg.solve(jac,b)

print('Sensitivities: ', x)

# Update state variables and loads of the system - assumed step size of 1.0
step = 1.0

teta1 = teta1 + step*x[0]
teta2 = teta2 + step*x[1]
v1 = v1 + step*x[2]
v2 = v2 + step*x[3]

P01 = P01 + step*beta[0]  # Presumed S = 1
P02 = P02 + step*beta[1]
Q01 = Q01 + step*alfa[0]
Q02 = Q02 + step*alfa[1]

if contpar == 'load':

    i = 0
    #----------------------------------------------------------------- Correction iterations  -------------------------------------
    #  Iterative loop for corrections
    while i < 3:
        # Calculate net injections
        print('Correction iteration: ', i + 1)
        p1, p2, p3, q1, q2, q3 = PQ_injections(v1, v2, v3, teta1, teta2, teta3)
        print(' ')
        #   print('Injction: ', pA1,p2,p3,q1,q2,q3)
        print('Net injection   P1 :', '{:9.5f}'.format(p1), '  P2 :','{:9.5f}'.format(p2), '  P3 :', '{:9.5f}'.format(p3), '\n',
              '               Q1 :', '{:9.5f}'.format(q1),  '  Q2 :', '{:9.5f}'.format(q2),  '  Q3 :', '{:9.5f}'.format(q3))
        f.write('\n' + 'Net injection   P1 :'+ '{:9.5f}'.format(p1) + '  P2 :'+'{:9.5f}'.format(p2) + '  P3 :'+ '{:9.5f}'.format(p3)+ '\n'+
                '                Q1 :'+ '{:9.5f}'.format(q1) +  '  Q2 :'+ '{:9.5f}'.format(q2)+  '  Q3 :'+ '{:9.5f}'.format(q3) +'\n\n')
        print(' ')

        #Jacobi elements
        #Active power
        dp1dt1, dp1dt2, dp1dt3, dp2dt1, dp2dt2, dp2dt3, dp3dt1, dp3dt2, dp3dt3 = dPdTeta(v1, v2, v3, teta1, teta2, teta3)
        dp1dv1, dp1dv2, dp1dv3, dp2dv1, dp2dv2, dp2dv3, dp3dv1, dp3dv2, dp3dv3 = dPdV(v1, v2, v3, teta1, teta2, teta3)

        # Reactive power
        dq1dt1, dq1dt2, dq1dt3, dq2dt1, dq2dt2, dq2dt3, dq3dt1, dq3dt2, dq3dt3 = dQdTeta(v1, v2, v3, teta1, teta2, teta3)
        dq1dv1, dq1dv2, dq1dv3, dq2dv1, dq2dv2, dq2dv3, dq3dv1, dq3dv2, dq3dv3 = dQdV(v1, v2, v3, teta1, teta2, teta3)



        # Display the Jacobian elements (screen and file
        #    print_jacobian_to_screen()
        #   write_jacobian_to_file()

        jac = np.array(
            [[dp1dt1, dp1dt2, dp1dv1, dp1dv2, beta[0]], [dp2dt1, dp2dt2, dp2dv1, dp2dv2, beta[1]], [dq1dt1, dq1dt2, dq1dv1, dq1dv2, alfa[0]],
             [dq2dt1, dq2dt2, dq2dv1, dq2dv2, alfa[1]], [0.0, 0.0, 0.0, 0.0, 1.0]])
        b = np.array([-p1-P01, -p2-P02, -q1-Q01, -q2-Q02, 0])
        #    print(' ')
        #    print('b-vector: ', b)
        print('RHS-vector   dp1 :', '{:10.7f}'.format(-p1-P01), '  dp2 :','{:10.7f}'.format(-p2-P02), '\n',
              '            dq1 :', '{:10.7f}'.format(-q1-Q01),  '  dq2 :', '{:10.7f}'.format(-q2-Q02))
        f.write('\n' +'RHS-vector   dp1 :'+ '{:10.7f}'.format(-p1-P01)+ '  dp2 :'+'{:10.7f}'.format(-p2-P02)+ '\n'+
                '             dq1 :'+ '{:10.7f}'.format(-q1-Q01)+  '  dq2 :'+ '{:10.7f}'.format(-q2-Q02) + '\n\n')
        print(' ')

        x = np.linalg.solve(jac,b)

        #   print(' ')
        #    print('x-vector: ', x)
        print('Correction   dt1 :', '{:10.7f}'.format(x[0]), '  dt2 :','{:10.7f}'.format(x[1]), '\n',


              '            dv1 :', '{:10.7f}'.format(x[2]),  '  dv2 :', '{:10.7f}'.format(x[3]))
        f.write('\n' + 'Correction   dt1 :'+ '{:10.7f}'.format(x[0])+ '  dt2 :'+'{:10.7f}'.format(x[1])+ '\n'+
                '             dv1 :'+ '{:10.7f}'.format(x[2])+  '  dv2 :'+ '{:10.7f}'.format(x[3]) + '\n\n')
        print(' ')

        # Update the voltage angles and magnitudes
        teta1 = teta1 + x[0]
        teta2 = teta2 + x[1]
        v1 = v1 + x[2]
        v2 = v2 + x[3]
        i = i + 1
elif contpar == 'volt':

    i = 0
    #----------------------------------------------------------------- Correction iterations  -------------------------------------
    #  Iterative loop for corrections
    while i < 3:
        # Calculate net injections
        print('Correction iteration: ', i + 1)
        p1, p2, p3, q1, q2, q3 = PQ_injections(v1, v2, v3, teta1, teta2, teta3)
        print(' ')
        #   print('Injction: ', pA1,p2,p3,q1,q2,q3)
        print('Net injection   P1 :', '{:9.5f}'.format(p1), '  P2 :','{:9.5f}'.format(p2), '  P3 :', '{:9.5f}'.format(p3), '\n',
              '               Q1 :', '{:9.5f}'.format(q1),  '  Q2 :', '{:9.5f}'.format(q2),  '  Q3 :', '{:9.5f}'.format(q3))
        f.write('\n' + 'Net injection   P1 :'+ '{:9.5f}'.format(p1) + '  P2 :'+'{:9.5f}'.format(p2) + '  P3 :'+ '{:9.5f}'.format(p3)+ '\n'+
                '                Q1 :'+ '{:9.5f}'.format(q1) +  '  Q2 :'+ '{:9.5f}'.format(q2)+  '  Q3 :'+ '{:9.5f}'.format(q3) +'\n\n')
        print(' ')

        #Jacobi elements
        #Active power
        dp1dt1, dp1dt2, dp1dt3, dp2dt1, dp2dt2, dp2dt3, dp3dt1, dp3dt2, dp3dt3 = dPdTeta(v1, v2, v3, teta1, teta2, teta3)
        dp1dv1, dp1dv2, dp1dv3, dp2dv1, dp2dv2, dp2dv3, dp3dv1, dp3dv2, dp3dv3 = dPdV(v1, v2, v3, teta1, teta2, teta3)

        # Reactive power
        dq1dt1, dq1dt2, dq1dt3, dq2dt1, dq2dt2, dq2dt3, dq3dt1, dq3dt2, dq3dt3 = dQdTeta(v1, v2, v3, teta1, teta2, teta3)
        dq1dv1, dq1dv2, dq1dv3, dq2dv1, dq2dv2, dq2dv3, dq3dv1, dq3dv2, dq3dv3 = dQdV(v1, v2, v3, teta1, teta2, teta3)



        # Display the Jacobian elements (screen and file
        #    print_jacobian_to_screen()
        #   write_jacobian_to_file()

        jac = np.array(
            [[dp1dt1, dp1dt2, dp1dv1, dp1dv2, beta[0]], [dp2dt1, dp2dt2, dp2dv1, dp2dv2, beta[1]], [dq1dt1, dq1dt2, dq1dv1, dq1dv2, alfa[0]],
             [dq2dt1, dq2dt2, dq2dv1, dq2dv2, alfa[1]], [0.0, 0.0, 0.0, 1.0, 0.0]])
        b = np.array([-p1-P01, -p2-P02, -q1-Q01, -q2-Q02, 0])
        #    print(' ')
        #    print('b-vector: ', b)
        print('RHS-vector   dp1 :', '{:10.7f}'.format(-p1-P01), '  dp2 :','{:10.7f}'.format(-p2-P02), '\n',
              '            dq1 :', '{:10.7f}'.format(-q1-Q01),  '  dq2 :', '{:10.7f}'.format(-q2-Q02))
        f.write('\n' +'RHS-vector   dp1 :'+ '{:10.7f}'.format(-p1-P01)+ '  dp2 :'+'{:10.7f}'.format(-p2-P02)+ '\n'+
                '             dq1 :'+ '{:10.7f}'.format(-q1-Q01)+  '  dq2 :'+ '{:10.7f}'.format(-q2-Q02) + '\n\n')
        print(' ')

        x = np.linalg.solve(jac,b)

        #   print(' ')
        #    print('x-vector: ', x)
        print('Correction   dt1 :', '{:10.7f}'.format(x[0]), '  dt2 :','{:10.7f}'.format(x[1]), '\n',


              '            dv1 :', '{:10.7f}'.format(x[2]),  '  dv2 :', '{:10.7f}'.format(x[3]))
        f.write('\n' + 'Correction   dt1 :'+ '{:10.7f}'.format(x[0])+ '  dt2 :'+'{:10.7f}'.format(x[1])+ '\n'+
                '             dv1 :'+ '{:10.7f}'.format(x[2])+  '  dv2 :'+ '{:10.7f}'.format(x[3]) + '\n\n')
        print(' ')

        # Update the voltage angles and magnitudes
        teta1 = teta1 + x[0]
        teta2 = teta2 + x[1]
        v1 = v1 + x[2]
        v2 = v2 + x[3]
        # Correct the load
        S = x[4]
        P01 = P01 + step*beta[0]*S  # Presumed S take the needed value
        P02 = P02 + step*beta[1]*S
        Q01 = Q01 + step*alfa[0]*S
        Q02 = Q02 + step*alfa[1]*S
        i = i + 1


p1, p2, p3, q1, q2, q3 = PQ_injections(v1, v2, v3, teta1, teta2, teta3)
print('Final solution: ','\n')

print('Net injection   P1 :', '{:9.5f}'.format(p1), '  P2 :','{:9.5f}'.format(p2), '  P3 :', '{:9.5f}'.format(p3), '\n',
      '               Q1 :', '{:9.5f}'.format(q1),  '  Q2 :', '{:9.5f}'.format(q2),  '  Q3 :', '{:9.5f}'.format(q3))
f.write('\n' + 'Net injection   P1 :'+ '{:9.5f}'.format(p1) + '  P2 :'+'{:9.5f}'.format(p2) + '  P3 :'+ '{:9.5f}'.format(p3)+ '\n'+
        '                Q1 :'+ '{:9.5f}'.format(q1) +  '  Q2 :'+ '{:9.5f}'.format(q2)+  '  Q3 :'+ '{:9.5f}'.format(q3) +'\n\n')
print(' ')

# Close the file
f.close()
