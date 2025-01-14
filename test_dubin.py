# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import leastsq
# import casadi

import math

# print('LSR')

#######################Dubins:LSR
# p1 = np.array([10, 10, 0*np.pi])
# p2 = np.array([25, 25, -0*np.pi])

#######################Dubins:RSR
# p1 = np.array([0, 0, 90*np.pi/180])
# p2 = np.array([15, 15, 0*np.pi/180])

#######################Dubins:RLR
# p1 = np.array([10, 10, 0*np.pi/180])
# p2 = np.array([15, 15, 180*np.pi/180])

#######################Dubins:LRL
# p1 = np.array([10, 10, 180*np.pi/180])
# p2 = np.array([15, 15, 0*np.pi/180])



# dp=p2-p1
# r=5
# d=np.linalg.norm(dp[0:2])/r


# #坐标变换
# theta=np.mod(np.arctan2(dp[1],dp[0]),2*np.pi)
# alpha=np.mod(p1[2]-theta,2*np.pi)
# beta=np.mod(p2[2]-theta,2*np.pi)

# TYPE=['LSL','LSR','RSL','RSR','RLR','LRL']

#######LSL

def LSL(alpha,beta,d):
    tmp0 = d + np.sin(alpha) - np.sin(beta)
    p_squared = 2 + (d*d) -(2*np.cos(alpha - beta)) + (2*d*(np.sin(alpha) - np.sin(beta)))
    if p_squared < 0:
        L = [np.inf, np.inf, np.inf, np.inf,]
    else:
        tmp1 =np.arctan2( (np.cos(beta)-np.cos(alpha)), tmp0)
        t = np.mod((-alpha + tmp1 ), 2*np.pi)
        p = np.sqrt(p_squared)
        q = np.mod((beta - tmp1 ), 2*np.pi)
        L=[t+p+q, t, p, q]
    return L

#######LSR
def LSR(alpha,beta,d):
        # tmp0 =  d + np.sin( alpha) - np.sin( beta)
        p_squared = -2 + ( d* d) + (2*np.cos( alpha -  beta)) + (2* d*(np.sin( alpha)+np.sin( beta)))
        # p_squared = 2 + ( d* d) -(2*np.cos( alpha -  beta)) + (2* d*(np.sin( alpha) - np.sin( beta)))
        if p_squared < 0:
            L = [np.inf, np.inf, np.inf, np.inf,]
        else:
            p = np.sqrt(p_squared)
            tmp2 = np.arctan2( (-np.cos( alpha)-np.cos( beta)), ( d+np.sin( alpha)+np.sin( beta)) ) - np.arctan2(-2.0, p)

            t = np.mod((- alpha + tmp2 ), 2*np.pi)
            
            q = np.mod((-np.mod(( beta), 2*np.pi) + tmp2 ), 2*np.pi)
            L=[t+p+q, t, p, q]
        return L


#######RSL
def RSL(alpha,beta,d):
        # tmp0 =  d + np.sin( alpha) - np.sin( beta)
        p_squared = (d*d) -2 + (2*np.cos(alpha - beta)) - (2*d*(np.sin(alpha)+np.sin(beta)))

        if p_squared < 0:
            L = [np.inf, np.inf, np.inf, np.inf,]
        else:
            p = np.sqrt(p_squared)
            tmp2 = np.arctan2( (np.cos(alpha)+np.cos(beta)), (d-np.sin(alpha)-np.sin(beta)) ) - np.arctan2(2.0, p)

            t = np.mod((alpha - tmp2), 2*np.pi)
            
            q = np.mod((beta - tmp2), 2*np.pi)
            L=[t+p+q, t, p, q]
        return L

#######RSR
def RSR(alpha,beta,d):
        tmp0 = d-np.sin(alpha)+np.sin(beta)
        p_squared = 2 + (d*d) -(2*np.cos(alpha - beta)) + (2*d*(np.sin(beta)-np.sin(alpha)))

        if p_squared < 0:
            L = [np.inf, np.inf, np.inf, np.inf,]
        else:
            tmp1 = np.arctan2( (np.cos(alpha)-np.cos(beta)), tmp0 )
            t=np.mod(( alpha - tmp1 ), 2*np.pi)
            p = np.sqrt(p_squared)
            q = np.mod((-beta + tmp1), 2*np.pi)
            L=[t+p+q, t, p, q]
        return L

#######RLR
def RLR(alpha,beta,d):
        tmp_rlr = (6.0 - d*d + 2*np.cos(alpha - beta) + 2*d*(np.sin(alpha)-np.sin(beta))) / 8.0
        # p_squared = 2 + (d*d) -(2*np.cos(alpha - beta)) + (2*d*(np.sin(beta)-np.sin(alpha)))

        if abs(tmp_rlr)>1 :
            L = [np.inf, np.inf, np.inf, np.inf,]
        else:
            p=np.mod(( 2*np.pi - np.arccos( tmp_rlr ) ), 2*np.pi)
            t = np.mod((alpha - np.arctan2( np.cos(alpha)- np.cos(beta), d- np.sin(alpha)+np.sin(beta) ) + np.mod(p/2, 2*np.pi)), 2*np.pi)
            q = np.mod((alpha - beta - t + np.mod(p, 2*np.pi)), 2*np.pi)
            L=[t+p+q, t, p, q]
        return L


#######RLR
def LRL(alpha,beta,d):
        tmp_rlr = (6.0 - d*d + 2*np.cos(alpha - beta) + 2*d*(-np.sin(alpha)+np.sin(beta))) / 8.0
        # p_squared = 2 + (d*d) -(2*np.cos(alpha - beta)) + (2*d*(np.sin(beta)-np.sin(alpha)))

        if abs(tmp_rlr)>1 :
            L = [np.inf, np.inf, np.inf, np.inf,]
        else:
            p=np.mod(( 2*np.pi - np.arccos( tmp_rlr ) ), 2*np.pi)
            t = np.mod((-alpha - np.arctan2( np.cos(alpha)- np.cos(beta), d+ np.sin(alpha)-np.sin(beta) ) + p/2), 2*np.pi)
            q = np.mod((np.mod(beta, 2*np.pi) - alpha -t + np.mod(p, 2*np.pi)), 2*np.pi)
            L=[t+p+q, t, p, q]
        return L

#计算每一段的终点, seg_param: 前进的单位长度，
def dubins_segment(seg_param, seg_init, seg_type):
    seg_end=[0]*3
    if seg_type=='L':
        seg_end[0] = seg_init[0] + np.sin(seg_init[2]+seg_param) - np.sin(seg_init[2])
        seg_end[1] = seg_init[1] - np.cos(seg_init[2]+seg_param) + np.cos(seg_init[2])
        seg_end[2] = seg_init[2] + seg_param

    elif seg_type == 'R':
        seg_end[0] = seg_init[0] - np.sin(seg_init[2]-seg_param) + np.sin(seg_init[2])
        seg_end[1] = seg_init[1] + np.cos(seg_init[2]-seg_param) - np.cos(seg_init[2])
        seg_end[2] = seg_init[2] - seg_param

    elif seg_type == 'S':
        seg_end[0] = seg_init[0] + np.cos(seg_init[2]) * seg_param
        seg_end[1] = seg_init[1] + np.sin(seg_init[2]) * seg_param
        seg_end[2] = seg_init[2]
    return seg_end



# L=LSR(alpha,beta,d)
# L=RSL(alpha,beta,d)
# L=RSR(alpha,beta,d)
# L=RLR(alpha,beta,d)
# L=LRL(alpha,beta,d)



# p_start = [0,  0,  p1[2]]
# # 第一段的终点
# mid1 = dubins_segment(L[1], p_start,'L')
# #end point of the second part
# mid2 = dubins_segment(L[2], mid1,'R')


# path=[]
###############LSL
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     t = step / r
#     if t<L[1]:
#         end_pt = dubins_segment(t, p_start,'L')

#     elif t<L[1]+L[2]:
#         end_pt = dubins_segment( t-L[1],mid1,'S')
#     else:
#         end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'L')

#     end_pt[0] = end_pt[0] * r + p1[0]
#     end_pt[1] = end_pt[1] * r + p1[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     path.append(end_pt)

# path_array = np.array(path)

###############LSR
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     t = step / r
#     if t<L[1]:
#         end_pt = dubins_segment(t, p_start,'L')

#     elif t<L[1]+L[2]:
#         end_pt = dubins_segment( t-L[1],mid1,'S')
#     else:
#         end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'R')

#     end_pt[0] = end_pt[0] * r + p1[0]
#     end_pt[1] = end_pt[1] * r + p1[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     path.append(end_pt)


###############RSR
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     t = step / r
#     if t<L[1]:
#         end_pt = dubins_segment(t, p_start,'R')

#     elif t<L[1]+L[2]:
#         end_pt = dubins_segment( t-L[1],mid1,'S')
#     else:
#         end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'R')

#     end_pt[0] = end_pt[0] * r + p1[0]
#     end_pt[1] = end_pt[1] * r + p1[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     path.append(end_pt)

###############RLR
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     t = step / r
#     if t<L[1]:
#         end_pt = dubins_segment(t, p_start,'R')

#     elif t<L[1]+L[2]:
#         end_pt = dubins_segment( t-L[1],mid1,'L')
#     else:
#         end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'R')

#     end_pt[0] = end_pt[0] * r + p1[0]
#     end_pt[1] = end_pt[1] * r + p1[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     path.append(end_pt)


###############LRL
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     t = step / r
#     if t<L[1]:
#         end_pt = dubins_segment(t, p_start,'L')

#     elif t<L[1]+L[2]:
#         end_pt = dubins_segment( t-L[1],mid1,'R')
#     else:
#         end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'L')

#     end_pt[0] = end_pt[0] * r + p1[0]
#     end_pt[1] = end_pt[1] * r + p1[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     path.append(end_pt)

# path_array = np.array(path)

# path_array = np.array(path)
# plt.plot(p1[0],p1[1],'ro')
# plt.plot(p2[0],p2[1],'ro')
# plt.plot(path_array[:,0],path_array[:,1],'b')
# plt.show()
# print('pause')