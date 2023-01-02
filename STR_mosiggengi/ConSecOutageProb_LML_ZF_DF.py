import numpy as np
import random
import math

class Outage(con_out, sec_out):
    def __init__(self):
        super(Outage, self).__init__()
        self.iteration = 0
        self.R_T = 0
        self.R_E = 0
        self.rho_B = 0
        self.rho_E = 0
        self.Es_dBm = 0
        self.Noise_Var_B = 0
        self.Noise_Var_E = 0

    def avg_cal(self):
        R_S = self.R_T - self.R_E
        Es = 10^(self.Es_dBm/10)/1000
        sigma_E = math.sqrt(self.Noise_Var_E)
        sigma_B = math.sqrt(self.Noise_Var_B)
        avg_SNR_B = (2*Es) / (sigma_B)^2
        avg_SNR_E = (2*Es) / (sigma_E)^2
        return avg_SNR_B, avg_SNR_E

    def 