class Function:
    def __init__(self):
        pass
    def Eff_phy_quantity(a1,a2):
        if a1 == 0 or a2 == 0:
            raise ValueError("调和平均要求两个数都非零")
        Eff_mean = 1 / (1/a1 + 1/a2)
        return Eff_mean



