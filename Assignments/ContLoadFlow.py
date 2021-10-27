from LoadFlowLineDriven import *


class ContLoadFlow(LoadFlowLineDriven):

    def __init__(self, Buses, Lines):
        super(ContLoadFlow, self).__init__(Buses, Lines)

    # Continuation Load Flow Procedure
    def CPF(self):
        LoadFlowLineDriven.NR(self, itLimit=4, errorLimit=10**-5)


Buses, Lines = System_Setup()
CLF = ContLoadFlow(Buses, Lines)
CLF.CPF()
