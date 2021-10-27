import Assignment1
from Assignment1.LoadFlowLineDriven import LoadFlowLineDriven


class ContLoadFlow(LoadFlowLineDriven):

    def __init__(self, Buses, Lines):
        super(ContLoadFlow, self).__init__(Buses, Lines)


Buses, Lines = Assignment1.System_Setup
CLF = ContLoadFlow(Buses, Lines)
CLF.CPF()
