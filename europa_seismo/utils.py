import numpy as np

def Mw_to_Mo(Mw,units='Nm'):
    '''
    Converts moment magnitude to scalar moment

    positional arguments-------------------------------------
    Mw: The moment magnitude
    units (optional): units of the scalar moment. default is 'Nm'.
                      also can use 'dynecm'
    '''
    if units=='Nm':
        Mo = 10.0**(1.5*Mw + 9.1)
    elif units=='dynecm':
        Mo = 10.0**(1.5*Mw + 16.1)

    return Mo

def Mo_to_Mw(Mo,units='Nm'):
    '''
    Converts scalar moment to moment magnitude

    positional arguments-------------------------------------
    Mo: The scalar moment
    units (optional): units of the scalar moment. default is 'Nm'.
                      also can use 'dynecm'
    '''
    if units=='Nm':
        Mw = (2./3.)*(np.log10(Mo)-9.1)
    elif units=='dynecm':
        Mw = (2./3.)*(np.log10(Mo)-16.1)

    return Mw
