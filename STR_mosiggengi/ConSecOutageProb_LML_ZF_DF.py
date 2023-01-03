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

    def avg_var_cal(self):
        R_S = self.R_T - self.R_E
        Es = 10^(self.Es_dBm/10)/1000
        sigma_E = math.sqrt(self.Noise_Var_E)
        sigma_B = math.sqrt(self.Noise_Var_B)
        self.avg_SNR_B = (2*Es) / (sigma_B)^2
        self.avg_SNR_E = (2*Es) / (sigma_E)^2
        self.Noise_Var_E_v = 1 - self.rho_E^2
        self.Noise_Var_B_v = 1 - self.rho_B^2
        return self.avg_SNR_B, self.avg_SNR_E, self.Noise_Var_B_v, self.Noise_Var_E_v

    def TS_cha_ge(self):
        for k in range(self.iteration):
            H_B1 = (np.random.randn(2,1) + math.sqrt(-1)*np.random.randn(2,1) / math.sqrt(2))
            H_B2 = np.dot(self.rho_B, H_B1) + np.dot((np.random.randn(2,1) + math.sqrt(-1)*np.random.randn(2,1)), math.sqrt(self.Noise_Var_B_v))