import numpy as np
import math
import itertools
# State space
# RHO = [0, 1000]
# THETA = [-math.pi,math.pi]
# PHI   = [-math.pi,math.pi]
# OWNSTATE_X = [-200, 1000]
# OWNSTATE_Y = [-200, 1000] #[100.0,200.0,300.0,400.0]  #ft/s
# OWNSTATE_Sig=[-math.pi,math.pi]

# boundary=np.array([RHO, THETA, PHI, OWNSTATE_X,OWNSTATE_Y, OWNSTATE_Sig])
# number = [50, 10, 10, 50, 50, 10]

#Sampling time
DELTA=1

# Space of the system
RHO = [0, 200]
THETA = [-math.pi,math.pi]
PHI   = [-math.pi,math.pi]
OWNSTATE_X = [-100, 800]
OWNSTATE_Y = [-100, 800] #[100.0,200.0,300.0,400.0]  #ft/s
OWNSTATE_Sig=[-math.pi,math.pi]

boundary=np.array([RHO, THETA, PHI, OWNSTATE_X,OWNSTATE_Y, OWNSTATE_Sig])
# the number of regions in each dim
NUMBER = [50, 10, 10, 50, 50, 10]

# Input Space
U_O=[-0.2,0.2]
V_O=[6,9]
U_I=[-0.1,0.1]
V_I=[10]

boundary_u=np.array([U_O, V_O])
NUMBER_u = [3,3]


class AIRCRAFT():
    def __init__(self):
        self.state_variables = ['rho', 'theta', 'phi',
                                 'x_o', 'y_o', 'sigma_o']
        self.n = len(self.state_variables)
        # self.p = np.size(self.B,1)
        
    def set_spec(self):
        spec=1111
        return spec
    

class Input_spec():
    
    def __init__(self):
        self.partition = {}
        # self.specification = {}
        # self.goal = {}
        self.space = {}
        
        
        self.space['U_O']= U_O
        self.space['V_O']= V_O
        self.space['U_I']= U_I
        self.space['V_I']= V_I

        self.partition['boundary']=boundary_u
        self.partition['number'] = NUMBER_u
        self.partition['dim']=np.shape(boundary_u)[0]


class UAV_spec():
    
    def __init__(self):
        self.partition = {}
        # self.specification = {}
        self.goal = {}
        self.space = {}
        self.control = {}

        self.control['uMin'] = [-0.1, 6] # rad, velocity
        self.control['uMax'] = [0.1, 9]

        self.space['RHO']= RHO
        self.space['THETA']= THETA
        self.space['PHI']= PHI
        self.space['OWNSTATE_X']= OWNSTATE_X
        self.space['OWNSTATE_Y']= OWNSTATE_Y
        self.space['OWNSTATE_Sig']=OWNSTATE_Sig

        self.partition['boundary'] =np.array([RHO, THETA, PHI, OWNSTATE_X,OWNSTATE_Y, OWNSTATE_Sig])
        self.partition['number'] = NUMBER

        #target state
        self.goal['ownship_t'] = [np.array(['all', 'all', 'all', [410,430],[550,570],[math.pi/2 - 0.5, math.pi/2 + 0.5]], dtype='object')]
        
        # unsafe states
        self.critical = [np.array([[0, 2], 'all', 'all', 'all', 'all', 'all'], dtype='object')]