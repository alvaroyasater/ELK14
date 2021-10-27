from Assignment1.LoadFlowLineDriven import LoadFlowLineDriven
from Assignment1.System_Setup import System_Setup


class ContLoadFlow(LoadFlowLineDriven):

    def __init__(self, Buses, Lines):
        super(ContLoadFlow, self).__init__(Buses, Lines)

    # Continuation Load Flow Procedure
    def ContLoadFlow(self):
        LoadFlowLineDriven.NR(itLimit=4, errorLimit=10**-5)



Buses, Lines = System_Setup()
CLF = ContLoadFlow(Buses, Lines)
CLF.CPF()
