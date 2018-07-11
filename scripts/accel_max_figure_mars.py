import instaseis
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo import acceleration

#read instaseis database
db = instaseis.open_db('/Volumes/G-RAID with Thunderbolt/axisem_databases/Mars_EH45Tcold_5s_noatten_database')
print db

#Mw = 3 thrust at the equator
source = instaseis.Source(latitude=0.0,longitude=0.0,depth_in_m=300.0,m_rr=3.98e13,m_pp=-3.98e13)

#receiver at 60 degrees distance
receiver = instaseis.Receiver(latitude=0.0,longitude=60.0)

#get seismograms
st = db.get_seismograms(source,receiver,components='ZRT',kind='displacement')

st.plot()
