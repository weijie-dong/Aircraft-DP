import numpy as np
import math
import itertools

from abstraction import Abstraction,Abstraction_u
from partition_function import computeRegionIdx
import matplotlib.pyplot as plt
from dubins import Dubins_p

#parameters of the system 
own_s_0=[0,0,0]
own_s_T=[420,560,0.5*np.pi]
intru_s_0=[-120,90,0]
intru_s_T=[600,450,0*np.pi]

v_i=10
u_i=[-0.2,0.2]
r_i=v_i/u_i[1]

T=90
######## Generate the nominal point of intruder
nom_path_intru=Dubins_p(intru_s_0,intru_s_T,r_i)
indx1,L=nom_path_intru.dubins() # indx1: record the kind of dubins(LSL or others)
nom_path_intru.mid_point()
nom_points_intru=[]
for i in range(T+1):
    nom_points_intru.append(nom_path_intru.nex_pos_i(v_i*i))
nom_points_intru = np.array(nom_points_intru)


#transfer absolute to relative coordinate
def abso2rela(own_s, intru_s):
    own_s=np.array(own_s)
    intru_s=np.array(intru_s)
    temp_n=own_s[0:2]-intru_s[0:2]
    rho=np.linalg.norm(own_s[0:2]-intru_s[0:2])
    theta=np.arctan2(temp_n[1], temp_n[0])
    
    phi=own_s[2]-intru_s[2] #############这里必须要角度。。。。

    return np.array([rho,theta,phi,own_s[0],own_s[1],own_s[2]])



# test=abso2rela(own_s_0, intru_s_0)
# print(test)
# print('Hello')
        
Ab_s=Abstraction()
# Define states of abstraction 划分系统状态空间
Ab_s.partion_state()
# # Initialize results dictionaries
# Ab_s.initialize_results()
# Define actions

Ab_u=Abstraction_u()
Ab_u.partion_state()
Ab_u.initialize_results()
# print(Ab_u.spec.partition['width'])

idx_all={(0, 0, 0, 0, 0, 0):0, (0, 0, 0, 0, 1, 0):1}



###dynamic of intruder
def dynaminc_intruder(state, input_i):
    '''
    state：intruder current state, np.array

    output：control input under Dubins path
    '''

    x_i_n=state[0]+input_i[1]*np.cos(state[2])
    y_i_n=state[1]+input_i[1]*np.sin(state[2])
    sigma_i_n=state[2]+input_i[0]

    pos_i_n=np.array([x_i_n,y_i_n,sigma_i_n])
    return pos_i_n

def dynamic_own(state,input_o):
    '''
    state: np.array()，只有own的状态值
    '''
    #next pos of ownship
    # pos=state[3:6]
    x_o_n=state[0]+input_o[1]*np.cos(state[2])
    y_o_n=state[1]+input_o[1]*np.sin(state[2])
    sigma_o_n=state[2]+input_o[0]

    pos_o_n=np.array([x_o_n,y_o_n,sigma_o_n])
    return pos_o_n

def dynamic(state,input_o):
    '''
    state: np.array()
    '''
    #intruder 的input，后续需改为dubins path的结果
    input_i=[[u_i[0],v_i],[u_i[1],v_i]]
    r_i=v_i/u_i[1]

    rho=state[0]
    theta=state[1]
    phi=state[2]
    pos=state[3:6]

    #current state of intruder
    pos_i=np.array([pos[0]+rho*np.cos(theta), pos[1]+rho*np.sin(theta), phi+pos[2]])

    #next pos of ownship
    x_o_n=pos[0]+input_o[1]*np.cos(pos[2])
    y_o_n=pos[1]+input_o[1]*np.sin(pos[2])
    sigma_o_n=pos[2]+input_o[0]

    pos_o_n=np.array([x_o_n,y_o_n,sigma_o_n])

    #next pos of intruder ship#######################
    # x_i_n=pos_i[0]+input_i[1]*np.cos(pos_i[2])
    # y_i_n=pos_i[1]+input_i[1]*np.sin(pos_i[2])
    # sigma_i_n=pos_i[2]+input_i[0]

    # pos_i_n=

    ############# three types of input
    pos_i_n=dict()
    #[lower, up,dubins]
    #two extrem input##############
    for i,inp in enumerate(input_i):
        pos_i_n[i]=dynaminc_intruder(pos_i, inp)

    ####Dubins point
    pT=np.array([600, 450, 0*np.pi/180])
    dub=Dubins_p(pos_i,pT,r_i)
    indx1,L=dub.dubins()
    dub.mid_point()
    nx_p=dub.nex_pos_i(v_i)
    pos_i_n[2]=np.array(nx_p)

    ######## three global state
    pos_n=[]
    for i in range(3):
        temp_n=pos_o_n[0:2]-pos_i_n[i][0:2]

        rho_n=np.linalg.norm(pos_o_n[0:2]-pos_i_n[i][0:2])
        theta_n=np.arctan2(temp_n[1], temp_n[0])
        u_i_delta=pos_i_n[i][2]-pos_i[2]
        phi_n=phi+u_i_delta-input_o[0] #############这里必须要角度。。。。

        pos_n.append(np.array([rho_n,theta_n,phi_n,x_o_n,y_o_n,sigma_o_n]))

    return pos_n

# compute transition probability
state_transition_dict = dict()

def transition_pro(Ab_s,Ab_u):
    print('begin compute tran-pro... \n')
    state_transition_dict = dict()
    for slice, index_t in Ab_s.partition['R']['idx'].items():
        state_transition_dict[index_t]=dict()
        pos_c=Ab_s.partition['R']['center'][index_t]
        for slice_u, index_t_u in Ab_u.partition['R']['idx'].items():
            state_transition_dict[index_t][index_t_u]=dict()
            input=Ab_u.partition['R']['center'][index_t_u]
            pos_n=np.array(dynamic(pos_c,input)) # 3*6

            #求下个点所在的状态序号，输出还是每一个维度的划分序号的组合
            indices_nonneg=[]
            for p in pos_n:
                _, indices_nonneg_t =computeRegionIdx(p,Ab_s.spec.partition, True)
                indices_nonneg.append(indices_nonneg_t)
                # print('\n indices_nonneg:',indices_nonneg)

            # transfer to tuple
            index_t_n=[]
            for v in indices_nonneg:
                # m.append(tuple(v)) 
                index_t_n.append(Ab_s.partition['R']['idx'][tuple(v)])
                
            # print('\n index_t_n:',index_t_n)
            
            # idx=str(m)
            PROB=[0.2,0.2,0.6]
            for i,sta in enumerate(index_t_n):
                if sta in state_transition_dict[index_t][index_t_u]:
                    state_transition_dict[index_t][index_t_u][sta] += PROB[i]
                else:
                    state_transition_dict[index_t][index_t_u][sta] = PROB[i]
            # print(state_transition_dict)
    
    return state_transition_dict
########################################################################
#计算转移概率
state_transition_dict=transition_pro(Ab_s,Ab_u)
#######################################################################

################### label function
def indicate_f(region_lis, point):
    if point in region_lis:
        return True #存在
    else:
        return False

################### 值迭代过程
# l=len(Ab_s.partitionp['R']['center'])
# 初始化状态值函数
# V = {s: 0 for s in Ab_s.partition['R']['idx_inv'].keys()}
print('\n ************** initialize value table')

l=len(Ab_s.partition['R']['idx_inv'])
V=np.zeros((T+1, l))
print(V.shape)
policy = np.zeros((T+1, l), dtype=int)

#initialize the value function
# assign the value table as 1 at goal state
for i, v in enumerate(V[T,:]):
    region_lis=Ab_s.partition['goal']
    if indicate_f(region_lis, i):
        V[T][i]=1
        # print('\n label:', i)


goal=Ab_s.partition['goal']
unsafe=Ab_s.partition['critical']

print('\n ************** begin value iteration')
for t in range(T-1,-1,-1):
    for s in Ab_s.partition['R']['idx_inv'].keys():
        Q = {}
        for a in Ab_u.partition['R']['idx_inv'].keys():
            sum=0
            for s2 in state_transition_dict[s][a].keys():
                # print('\n s_2:', s2)
                # print('\n V[T][s2]:', V[T][s2])
                # print('\n state_transition_dict[s][a][s2]:',state_transition_dict[s][a][s2])
                sum +=V[t+1][s2]*state_transition_dict[s][a][s2]
            Q[a]=indicate_f(goal, s) + (not indicate_f(goal, s)) * (not indicate_f(unsafe, s))*sum
            # if Q[a] != 0:
            #     print('not zero')
        temp=max(Q.values())
        V[t,s]=temp
        max_keys = [key for key, value in Q.items() if value == temp]
        policy[t, s] = max_keys[0]
        # if temp != 0:
        #         print('not zero')
        # print('Hello')
###############################


########################################## draw the online control trajactory
Path_own=[own_s_0]
s_own_c=own_s_0
for i in range(T+1):
    s_abs_c=abso2rela(s_own_c, nom_points_intru[i]) #绝对坐标转换后的相对坐标值
    _, indices_nonneg_t =computeRegionIdx(s_abs_c,Ab_s.spec.partition, True)
    s_rela_c=Ab_s.partition['R']['idx'][tuple(indices_nonneg_t)] #相对坐标在partition中的绝对序号
    # 获取当前的policy: np.zeros((T, l), dtype=int)
    u_abs_c=policy[i,s_rela_c] # input 的绝对序号值
    u_c=Ab_u.partition['R']['center'][u_abs_c]# input 的实际值
    s_own_c=dynamic_own(s_abs_c,u_c)
    Path_own.append(s_own_c)


traj_own = np.array(Path_own)

plt.plot(own_s_0[0],own_s_0[1],'go')
plt.plot(own_s_T[0],own_s_T[1],'go')
plt.plot(intru_s_T[0],intru_s_T[1],'ro')
plt.plot(intru_s_0[0],intru_s_0[1],'ro')
plt.plot(nom_points_intru[:,0],nom_points_intru[:,1],'x')
plt.plot(traj_own[:,0],traj_own[:,1],'-')
plt.show()
print('Hello')
