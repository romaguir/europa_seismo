import numpy as np
import scipy.stats as st

class gr_distribution(st.rv_continuous):
    def _pdf(self,b,M,Ntot):
        return Ntot*10**(-b*M)
        
