import europa_seismo
from europa_seismo.europa_seismo import minos

#runmodel
modesfile='modes.out'
eigfile='eig.out'
eps=1e-10
wgrav=1.315
jcom=3
lmin=0
lmax=1000
wmin=0.1
wmax=100.0
nmin=0.0
nmax=0.0
modelspath = europa_seismo.__path__[0] +'/data/models/'
modelfile = modelspath + 'aicehot_10km_prem_SIMPLE.txt'
print modelfile

#run model
minos.main(modelfile,modesfile,eigfile,eps,wgrav,jcom,lmin,lmax,wmin,wmax,nmin,nmax)
