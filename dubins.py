# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import leastsq
# import casadi

import math

from test_dubin import LSL, LSR, RSL, RSR, RLR, LRL

'''
###############################################################################
#生成圆周围的点
def get_circle_pnt(u, pos, radius, offset=0.0):
    x = pos[0] + radius * np.cos(2 * np.pi * u + offset)
    y = pos[1] + radius * np.sin(2 * np.pi * u + offset)
    return x, y

def get_circle(pos, radius, offset=0.0):
    u = np.linspace(0, 1, 50)#产生50个点，即将2 pi划分成50份
    return get_circle_pnt(u, pos, radius, offset)

#获得旋转矩阵，把原始系统角度变为0度
def get_frame(angle):
    return np.array([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])

radius = 5
wp_1 = np.array([2, 3])
heading_1 = 10 * np.pi / 180 #heading表示朝向的角度
wp_2 = np.array([20, 32])
heading_2 = 130 * np.pi / 180

frame_1 = get_frame(heading_1)
frame_2 = get_frame(heading_2)

m=frame_1[:, 1]
#计算相切的两个圆
center_1 = dict(R=wp_1 - radius * frame_1[:, 1].flatten(),
                L=wp_1 + radius * frame_1[:, 1].flatten())

center_2 = dict(R=wp_2 - radius * frame_2[:, 1].flatten(),
                L=wp_2 + radius * frame_2[:, 1].flatten())

def plot_base():
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)

    # Plot first waypoint's L and R circles
    #沿着前进方向画长度为5的线
    ax.plot([wp_1[0], wp_1[0] + 5 * np.cos(heading_1)], [wp_1[1], wp_1[1] + 5 * np.sin(heading_1)], color='xkcd:neon pink', linewidth=3)
    ax.plot([wp_1[0]], [wp_1[1]], '.r', markersize=10)

    #围绕圆心画圆
    x, y = get_circle(center_1['L'], radius)
    ax.plot(x, y, '-.r')
    x, y = get_circle(center_1['R'], radius)
    ax.plot(x, y, '--r')

    # Plot second waypoint's L and R circles
    x, y = get_circle(wp_2, radius)
    ax.plot([wp_2[0], wp_2[0] + 5 * np.cos(heading_2)], [wp_2[1], wp_2[1] + 5 * np.sin(heading_2)], color='xkcd:neon pink', linewidth=3)
    ax.plot([wp_2[0]], [wp_2[1]], '.b', markersize=10)

    x, y = get_circle(center_2['L'], radius)
    ax.plot(x, y, '-.b')
    x, y = get_circle(center_2['R'], radius)
    ax.plot(x, y, '--b')

    ax.axis('equal')
    ax.grid(True)
    fig.show()
    return ax

# f=plot_base()
# f.show()
print('pause')
#############################################################################
'''

# print('RSR')

# u1 = casadi.SX.sym('u1')
# u2 = casadi.SX.sym('u2')
# get_tangents(center_1['R'], radius, heading_1, delta_1=-1, center_2['R'], radius, heading_2, delta_2=-1)
# delta_1 的作用是什么？
'''
def get_tangents(center_1, radius_1, heading_1, delta_1, center_2, radius_2, heading_2, delta_2):
    output = dict()

    phi_1 = 2 * np.pi * u1 * delta_1 + heading_1 - delta_1 * np.pi / 2
    phi_2 = 2 * np.pi * u2 * delta_2 + heading_2 - delta_2 * np.pi / 2

    u1_func = lambda angle: (angle - heading_1 + delta_1 * np.pi / 2) / (delta_1 * 2 * np.pi)
    u2_func = lambda angle: (angle - heading_2 + delta_2 * np.pi / 2) / (delta_2 * 2 * np.pi)
    # Make tangents vector functions
    tan_1 = casadi.cross(np.array([0, 0, 1]), 
                         np.array([delta_1 * radius * np.cos(phi_1), delta_1 * radius * np.sin(phi_1), 0]))[0:2]
    tan_2 = casadi.cross(np.array([0, 0, 1]), 
                         np.array([delta_2 * radius * np.cos(phi_2), delta_2 * radius * np.sin(phi_2), 0]))[0:2]

    # Make circle functions
    circle_1_func = center_1 + radius_1 * np.array([np.cos(phi_1), np.sin(phi_1)])
    circle_2_func = center_2 + radius_2 * np.array([np.cos(phi_2), np.sin(phi_2)])

    # Plot the circles
    ax = plot_base()

    ## Plot the circle center points
    ax.plot([center_1[0]], [center_1[1]], marker='.', color='xkcd:baby blue', markersize=10)
    ax.plot([center_2[0]], [center_2[1]], marker='.', color='xkcd:baby blue', markersize=10)

    ## Plot a couple of tangent vectors
    for i in np.linspace(0.2, 0.8, 5):
        c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [i])
        t1 = casadi.substitute(tan_1, casadi.vertcat(*[u1]), [i])
        t1 *= 5 / casadi.norm_2(t1)
        ax.plot([c1[0], c1[0] + t1[0]],
                [c1[1], c1[1] + t1[1]], 'k')

        c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [i])
        t2 = casadi.substitute(tan_2, casadi.vertcat(*[u2]), [i])
        t2 *= 5 / casadi.norm_2(t2)
        ax.plot([c2[0], c2[0] + t2[0]],
                [c2[1], c2[1] + t2[1]], 'k')

    ## Plot line connecting the circle centers
    ax.plot([center_1[0], center_2[0]], [center_1[1], center_2[1]], linestyle='--', color='xkcd:lavender')

    ## Plot the starting point of each circle
    c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [0])
    c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [0])
    ax.plot([c1[0]], [c1[1]], marker='o', color='xkcd:reddish orange', markersize=15)
    ax.plot([c2[0]], [c2[1]], marker='o', color='xkcd:reddish orange', markersize=15)

    # Compute the line connecting the circle's centers
    d = center_2 - center_1
    # Calculate normal vector to the connecting line
    n = np.dot(get_frame(np.pi / 2), d / np.linalg.norm(d))

    ## Plotting the normal vectors
    ax.plot([center_1[0], center_1[0] + radius_1 * n[0]], [center_1[1], center_1[1] + radius_1 * n[1]], linestyle='--', color='xkcd:hot pink')
    ax.plot([center_2[0], center_2[0] + radius_2 * n[0]], [center_2[1], center_2[1] + radius_2 * n[1]], linestyle='--', color='xkcd:hot pink')

    ##########################################################
    # Compute the first tangent
    ## Compute the normal vector's angle
    n_angle = np.arctan2(n[1], n[0])    
    ## Compute the parameter for the tangent points on both circles
    u1_opt = u1_func(n_angle)    
    if u1_opt < 0:
        u1_opt = u1_opt + 1 
    u2_opt = u2_func(n_angle)
    if u2_opt < 0:
        u2_opt = u2_opt + 1

    ## Compute the points on the circles for the first tangent
    c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt])
    c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])

    tangent_1 = c2 - c1
    tangent_1 /= casadi.norm_2(tangent_1)

    ax.plot([c1[0], c2[0]], [c1[1], c2[1]], linestyle='--', color='xkcd:kelly green')

    ## Compute the tangent vectors on the circles
    t1 = casadi.substitute(tan_1, casadi.vertcat(*[u1]), [u1_opt])
    t1 /= casadi.norm_2(t1)
    t2 = casadi.substitute(tan_2, casadi.vertcat(*[u2]), [u2_opt])
    t2 /= casadi.norm_2(t2)

    diff = float(casadi.norm_2(tangent_1 - t1) + casadi.norm_2(tangent_1 - t2))

    if np.isclose(diff, 0):
        u = np.arange(0, u1_opt, 0.001)        
        output['C1'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [ui]) for ui in u]
        output['S'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt]), casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])]
        u = np.arange(u2_opt, 1, 0.001)
        output['C2'] = [casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [ui]) for ui in u]

    ## Plot the tangent vectors on the circles that are parallel to the first tangent 
    ax.plot([c1[0], c1[0] + radius_1 * t1[0]], [c1[1], c1[1] + radius_1 * t1[1]], linestyle='-', color='xkcd:bright purple')
    ax.plot([c2[0], c2[0] + radius_2 * t2[0]], [c2[1], c2[1] + radius_2 * t2[1]], linestyle='-', color='xkcd:bright purple')

    ##########################################################
    # Compute the second tangent    
    n_angle = np.arctan2(-n[1], -n[0])
    u1_opt = u1_func(n_angle)    
    if u1_opt < 0:
        u1_opt = u1_opt + 1 
    u2_opt = u2_func(n_angle)
    if u2_opt < 0:
        u2_opt = u2_opt + 1

    c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt])
    c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])

    tangent_2 = c2 - c1
    tangent_2 /= casadi.norm_2(tangent_2)

    ## Plotting the second tangent
    ax.plot([c1[0], c2[0]], [c1[1], c2[1]], linestyle='--', color='xkcd:kelly green')

    ## Compute the tangent vectors on the circles
    t1 = casadi.substitute(tan_1, casadi.vertcat(*[u1]), [u1_opt])
    t1 /= casadi.norm_2(t1)
    t2 = casadi.substitute(tan_2, casadi.vertcat(*[u2]), [u2_opt])
    t2 /= casadi.norm_2(t2)

    diff = float(casadi.norm_2(tangent_2 - t1) + casadi.norm_2(tangent_2 - t2))

    if np.isclose(diff, 0):
        u = np.arange(0, u1_opt, 0.001)        
        output['C1'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [ui]) for ui in u]
        output['S'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt]), casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])]
        u = np.arange(u2_opt, 1, 0.001)
        output['C2'] = [casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [ui]) for ui in u]

    ## Plot the tangent vectors on the circles that are parallel to the second tangent 
    ax.plot([c1[0], c1[0] + radius_1 * t1[0]], [c1[1], c1[1] + radius_1 * t1[1]], linestyle='-', color='xkcd:bright purple')
    ax.plot([c2[0], c2[0] + radius_2 * t2[0]], [c2[1], c2[1] + radius_2 * t2[1]], linestyle='-', color='xkcd:bright purple')

    ##########################################################
    # Computing inner tangents
    # Calculate the intersection point of the two tangent lines
    xp = (center_1[0] * radius_1 + center_2[0] * radius_2) / (radius_1 + radius_2)
    yp = (center_1[1] * radius_1 + center_2[1] * radius_2) / (radius_1 + radius_2)

    ax.plot([xp], [yp], '.r', markersize=10)
    # Third and fourth tangents
    xt1 = (radius_1**2 * (xp - center_1[0]) + radius_1 * (yp - center_1[1]) * np.sqrt((xp - center_1[0])**2 + (yp - center_1[1])**2 - radius_1**2)) / ((xp - center_1[0])**2 + (yp - center_1[1])**2) + center_1[0]
    xt2 = (radius_1**2 * (xp - center_1[0]) - radius_1 * (yp - center_1[1]) * np.sqrt((xp - center_1[0])**2 + (yp - center_1[1])**2 - radius_1**2)) / ((xp - center_1[0])**2 + (yp - center_1[1])**2) + center_1[0]

    yt1 = ((radius_1**2 * (yp - center_1[1])) - radius_1 * (xp - center_1[0]) * np.sqrt((xp - center_1[0])**2 + (yp - center_1[1])**2 - radius_1**2)) / ((xp - center_1[0])**2 + (yp - center_1[1])**2) + center_1[1]
    yt2 = ((radius_1**2 * (yp - center_1[1])) + radius_1 * (xp - center_1[0]) * np.sqrt((xp - center_1[0])**2 + (yp - center_1[1])**2 - radius_1**2)) / ((xp - center_1[0])**2 + (yp - center_1[1])**2) + center_1[1]

    ## Plotting the tangent points on the first circle
    ax.plot([xt1, xt2], [yt1, yt2], '.r', markersize=10)

    xt3 = (radius_2**2 * (xp - center_2[0]) + radius_2 * (yp - center_2[1]) * np.sqrt((xp - center_2[0])**2 + (yp - center_2[1])**2 - radius_2**2)) / ((xp - center_2[0])**2 + (yp - center_2[1])**2) + center_2[0]
    xt4 = (radius_2**2 * (xp - center_2[0]) - radius_2 * (yp - center_2[1]) * np.sqrt((xp - center_2[0])**2 + (yp - center_2[1])**2 - radius_2**2)) / ((xp - center_2[0])**2 + (yp - center_2[1])**2) + center_2[0]

    yt3 = ((radius_2**2 * (yp - center_2[1])) - radius_2 * (xp - center_2[0]) * np.sqrt((xp - center_2[0])**2 + (yp - center_2[1])**2 - radius_2**2)) / ((xp - center_2[0])**2 + (yp - center_2[1])**2) + center_2[1]
    yt4 = ((radius_2**2 * (yp - center_2[1])) + radius_2 * (xp - center_2[0]) * np.sqrt((xp - center_2[0])**2 + (yp - center_2[1])**2 - radius_2**2)) / ((xp - center_2[0])**2 + (yp - center_2[1])**2) + center_2[1]

    ## Plotting the tangent points on the second circle
    ax.plot([xt3, xt4], [yt3, yt4], '.r', markersize=10)

    # Third tangent
    u1_opt = u1_func(np.arctan2(yt1 - center_1[1], xt1 - center_1[0]))
    if u1_opt < 0:
        u1_opt = u1_opt + 1     
    u2_opt = u2_func(np.arctan2(yt3 - center_2[1], xt3 - center_2[0]))
    if u2_opt < 0:
        u2_opt = u2_opt + 1

    c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt])
    c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])

    t1 = casadi.substitute(tan_1, casadi.vertcat(*[u1]), [u1_opt])
    t1 /= casadi.norm_2(t1)
    t2 = casadi.substitute(tan_2, casadi.vertcat(*[u2]), [u2_opt])
    t2 /= casadi.norm_2(t2)

    ## Plot the tangent vectors on the circles that are parallel to the third tangent 
    ax.plot([c1[0], c1[0] + radius_1 * t1[0]], [c1[1], c1[1] + radius_1 * t1[1]], linestyle='-', color='xkcd:bright purple')
    ax.plot([c2[0], c2[0] + radius_2 * t2[0]], [c2[1], c2[1] + radius_2 * t2[1]], linestyle='-', color='xkcd:bright purple')

    tangent_3 = np.array([xt3 - xt1, yt3 - yt1])
    tangent_3 /= np.linalg.norm(tangent_3)

    diff = float(casadi.norm_2(tangent_3 - t1) + casadi.norm_2(tangent_3 - t2))

    if np.isclose(diff, 0):
        u = np.arange(0, u1_opt, 0.001)        
        output['C1'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [ui]) for ui in u]
        output['S'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt]), casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])]
        u = np.arange(u2_opt, 1, 0.001)
        output['C2'] = [casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [ui]) for ui in u]

    # Fourth tangent
    u1_opt = u1_func(np.arctan2(yt2 - center_1[1], xt2 - center_1[0]))
    if u1_opt < 0:
        u1_opt = u1_opt + 1
    u2_opt = u2_func(np.arctan2(yt4 - center_2[1], xt4 - center_2[0]))
    if u2_opt < 0:
        u2_opt = u2_opt + 1

    c1 = casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt])
    c2 = casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])

    t1 = casadi.substitute(tan_1, casadi.vertcat(*[u1]), [u1_opt])
    t1 /= casadi.norm_2(t1)
    t2 = casadi.substitute(tan_2, casadi.vertcat(*[u2]), [u2_opt])
    t2 /= casadi.norm_2(t2)

    ## Plot the tangent vectors on the circles that are parallel to the fourth tangent 
    ax.plot([c1[0], c1[0] + radius_1 * t1[0]], [c1[1], c1[1] + radius_1 * t1[1]], linestyle='-', color='xkcd:bright purple')
    ax.plot([c2[0], c2[0] + radius_2 * t2[0]], [c2[1], c2[1] + radius_2 * t2[1]], linestyle='-', color='xkcd:bright purple')

    tangent_4 = np.array([xt4 - xt2, yt4 - yt2])
    tangent_4 /= np.linalg.norm(tangent_4)

    diff = float(casadi.norm_2(tangent_4 - t1) + casadi.norm_2(tangent_4 - t2))

    if np.isclose(diff, 0):
        u = np.arange(0, u1_opt, 0.001)        
        output['C1'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [ui]) for ui in u]
        output['S'] = [casadi.substitute(circle_1_func, casadi.vertcat(*[u1]), [u1_opt]), casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [u2_opt])]
        u = np.arange(u2_opt, 1, 0.001)
        output['C2'] = [casadi.substitute(circle_2_func, casadi.vertcat(*[u2]), [ui]) for ui in u]

    ax.plot([xt1, xt3], [yt1, yt3], '--c')
    ax.plot([xt2, xt4], [yt2, yt4], '--c')

    #########################################################
    # Plot the path            
    ax.plot([x[0] for x in output['C1']], [x[1] for x in output['C1']], color='xkcd:golden yellow', linewidth=3)
    ax.plot([x[0] for x in output['S']], [x[1] for x in output['S']], color='xkcd:vermillion', linewidth=3)
    ax.plot([x[0] for x in output['C2']], [x[1] for x in output['C2']], color='xkcd:bright magenta', linewidth=3)

get_tangents(center_1['R'], radius, heading_1, -1, center_2['R'], radius, heading_2, -1)
'''

TYPE=['LSL','LSR','RSL','RSR','RLR','LRL']

class Dubins_p():
    
    def __init__(self,p0,pT,r=5):
        self.p0=np.array(p0)
        self.pT=np.array(pT)
        self.L=[]
        self.ind=0
        self.mid1=[]
        self.mid2=[]
        #坐标变换
        self.dp=self.pT-self.p0
        self.r=r
        self.d=np.linalg.norm(self.dp[0:2])/self.r
        self.theta=np.mod(np.arctan2(self.dp[1],self.dp[0]),2*np.pi)
        self.alpha=np.mod(p0[2]-self.theta,2*np.pi)
        self.beta=np.mod(pT[2]-self.theta,2*np.pi)

        self.p_start= [0,  0,  p0[2]]
        
    # TYPE=['LSL','LSR','RSL','RSR','RLR','LRL']
    def dubins(self):
        self.L.append(LSL(self.alpha,self.beta,self.d))
        self.L.append(LSR(self.alpha,self.beta,self.d))
        self.L.append(RSL(self.alpha,self.beta,self.d))
        self.L.append(RSR(self.alpha,self.beta,self.d))
        self.L.append(RLR(self.alpha,self.beta,self.d))
        self.L.append(LRL(self.alpha,self.beta,self.d))
        L_array=np.array(self.L)
        self.ind=np.argmin(L_array[:,0])
        return self.ind, L_array[self.ind,:]
    
    # def LSL(self):
    #     tmp0 = self.d + np.sin(self.alpha) - np.sin(self.beta)
    #     p_squared = 2 + (self.d*self.d) -(2*np.cos(self.alpha - self.beta)) + (2*self.d*(np.sin(self.alpha) - np.sin(self.beta)))
    #     if p_squared < 0:
    #         L = [np.inf, np.inf, np.inf, np.inf,]
    #     else:
    #         tmp1 =np.arctan2( (np.cos(self.beta)-np.cos(self.alpha)), tmp0)
    #         t = np.mod((-self.alpha + tmp1 ), 2*np.pi)
    #         p = np.sqrt(p_squared)
    #         q = np.mod((self.beta - tmp1 ), 2*np.pi)

    #     return [t+p+q, t, p, q]
    
    # def LSR(self):
    #     # tmp0 = self.d + np.sin(self.alpha) - np.sin(self.beta)
    #     p_squared = -2 + (self.d*self.d) + (2*np.cos(self.alpha - self.beta)) + (2*self.d*(np.sin(self.alpha)+np.sin(self.beta)))
    #     # p_squared = 2 + (self.d*self.d) -(2*np.cos(self.alpha - self.beta)) + (2*self.d*(np.sin(self.alpha) - np.sin(self.beta)))
    #     if p_squared < 0:
    #         L = [np.inf, np.inf, np.inf, np.inf,]
    #     else:
    #         p = np.sqrt(p_squared)
    #         tmp2 = np.arctan2( (-np.cos(self.alpha)-np.cos(self.beta)), (self.d+np.sin(self.alpha)+np.sin(self.beta)) ) - np.arctan2(-2.0, p)

    #         t = np.mod((-self.alpha + tmp2 ), 2*np.pi)
            
    #         q = np.mod((-np.mod((self.beta), 2*np.pi) + tmp2 ), 2*np.pi)

    #     return [t+p+q, t, p, q]
    
    def dubins_segment(self,seg_param, seg_init, seg_type):
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
    
    #计算两段距离的终点 
    # 要根据实际种类更改mid——point的内容
    def mid_point(self):
        #第一段的终点
        # a=self.L[self.ind][1]
        # b=self.p_start
        # c=TYPE[self.ind][0]
        # self.mid1 = self.dubins_segment(a, b)
        self.mid1 = self.dubins_segment(self.L[self.ind][1], self.p_start,TYPE[self.ind][0])
        #end point of the second part
        self.mid2 = self.dubins_segment(self.L[self.ind][2], self.mid1,TYPE[self.ind][1])


    def nex_pos_i(self,disp):
        '''
        init_p:list
        disp:float
        type:str
        '''
        t = disp / self.r
        L=self.L[self.ind]
        type=TYPE[self.ind]
        
        #########already arrived the target point
        if t>L[0]:
            end_pt = self.pT
            return end_pt

        start = [0,  0,  self.p0[2]]
        if t<L[1]:
            end_pt = self.dubins_segment(t, start,type[0])
        elif t<L[1]+L[2]:
            end_pt = self.dubins_segment(t-L[1],self.mid1,type[1])
        elif t<L[0]:
            end_pt = self.dubins_segment(t-(L[1]+L[2]),self.mid2,type[2])

        end_pt[0] = end_pt[0] * self.r + self.p0[0]
        end_pt[1] = end_pt[1] * self.r + self.p0[1]
        end_pt[2] = np.mod(end_pt[2], 2*np.pi)

        return end_pt

###dynamic of intruder
# def dynaminc_intruder(state, input_i):
#     '''
#     state：intruder current state, np.array

#     output：control input under Dubins path
#     '''

#     x_i_n=state[0]+input_i[1]*np.cos(state[2])
#     y_i_n=state[1]+input_i[1]*np.sin(state[2])
#     sigma_i_n=state[2]+input_i[0]

#     pos_i_n=np.array([x_i_n,y_i_n,sigma_i_n])
#     return pos_i_n

#######################Dubins:LSL
# p1 = np.array([10, 10, 0*np.pi])
# p2 = np.array([15, 15, -np.pi])

#######################Dubins:LSR
# p1 = np.array([-120, 90, 0*np.pi/180])
# p2 = np.array([600, 450, 0*np.pi/180])

##########upper input######################################

# #dubins##########################
# v=10
# u=0.1
# r=v/u
# dub=Dubins_p(p1,p2,r)
# indx1,L=dub.dubins()
# dub.mid_point()
# nx_p=dub.nex_pos_i(10)

# print("Dubins-type:", TYPE[dub.ind])
# # r=5
# path=[]
# for step in np.arange(0,L[0]*r+0.5,0.5):
#     e=dub.nex_pos_i(step)
#     path.append(e)

# path_array = np.array(path)

# path_array = np.array(path)
# plt.plot(p1[0],p1[1],'ro')
# plt.plot(p2[0],p2[1],'ro')
# plt.plot(path_array[:,0],path_array[:,1],'b')
# plt.plot(nx_p[0],nx_p[1],'gx')
# plt.show()
# print('pause')

##############################
# dp=p2-p1
# r=5
# d=np.linalg.norm(dp)/r


# #坐标变换
# theta=np.mod(np.arctan2(dp[1],dp[0]),2*np.pi)
# alpha=np.mod(p1[2]-theta,2*np.pi)
# beta=np.mod(p2[2]-theta,2*np.pi)

# TYPE=['LSL','LSR','RSL','RSR','RLR','LRL']
# L=dict()
# L['LSL']=

# def LSL(alpha,beta,d):
#     tmp0 = d + np.sin(alpha) - np.sin(beta)
#     p_squared = 2 + (d*d) -(2*np.cos(alpha - beta)) + (2*d*(np.sin(alpha) - np.sin(beta)))
#     if p_squared < 0:
#         L = [np.inf, np.inf, np.inf, np.inf,]
#     else:
#         tmp1 =np.arctan2( (np.cos(beta)-np.cos(alpha)), tmp0)
#         t = np.mod((-alpha + tmp1 ), 2*np.pi)
#         p = np.sqrt(p_squared)
#         q = np.mod((beta - tmp1 ), 2*np.pi)

#     return [t+p+q, t, p, q]
# # 选出最好的下一步的点
# L=LSL(alpha,beta,d)

# p_start = [0,  0,  p1[2]]

#计算每一段的终点, seg_param: 前进的单位长度，
# def dubins_segment(seg_param, seg_init, seg_type):
#     seg_end=[0]*3
#     if seg_type=='L':
#         seg_end[0] = seg_init[0] + np.sin(seg_init[2]+seg_param) - np.sin(seg_init[2])
#         seg_end[1] = seg_init[1] - np.cos(seg_init[2]+seg_param) + np.cos(seg_init[2])
#         seg_end[2] = seg_init[2] + seg_param

#     elif seg_type == 'R':
#         seg_end[0] = seg_init[0] - np.sin(seg_init[2]-seg_param) + np.sin(seg_init[2])
#         seg_end[1] = seg_init[1] + np.cos(seg_init[2]-seg_param) - np.cos(seg_init[2])
#         seg_end[2] = seg_init[2] - seg_param

#     elif seg_type == 'S':
#         seg_end[0] = seg_init[0] + np.cos(seg_init[2]) * seg_param
#         seg_end[1] = seg_init[1] + np.sin(seg_init[2]) * seg_param
#         seg_end[2] = seg_init[2]
#     return seg_end

#第一段的终点
# mid1 = dubins_segment(L[1], p_start,'L')
# #end point of the second part
# mid2 = dubins_segment(L[2], mid1,'S')

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


# def nex_pos_i(dub,disp,L,type,r):
#     '''
#     init_p:list
#     disp:float
#     type:str
#     r:float
#     '''
#     t = disp / r
#     start = [0,  0,  dub.p0[2]]
#     if t<L[1]:
#         end_pt = dub.dubins_segment(t, start,type[0])
#     elif t<L[1]+L[2]:
#         end_pt = dub.dubins_segment( t-L[1],dub.mid1,type[1])
#     else:
#         end_pt = dub.dubins_segment( t-(L[1]+L[2]),dub.mid2,type[2])

#     end_pt[0] = end_pt[0] * r + dub.p0[0]
#     end_pt[1] = end_pt[1] * r + dub.p0[1]
#     end_pt[2] = np.mod(end_pt[2], 2*np.pi)

#     return end_pt

'''
for step in np.arange(0,L[0]*r+0.5,0.5):
    e=nex_pos_i(p_start,step,L,'LSL',r)
    # t = step / r
    # if t<L[1]:
    #     end_pt = dubins_segment(t, p_start,'L')

    # elif t<L[1]+L[2]:
    #     end_pt = dubins_segment( t-L[1],mid1,'S')
    # else:
    #     end_pt = dubins_segment( t-(L[1]+L[2]),mid2,'L')

    # end_pt[0] = end_pt[0] * r + p1[0]
    # end_pt[1] = end_pt[1] * r + p1[1]
    # end_pt[2] = np.mod(end_pt[2], 2*np.pi)

    path.append(e)

path_array = np.array(path)


   

# print(path_array[:,0])
plt.plot(p1[0],p1[1],'ro')
plt.plot(p2[0],p2[1],'ro')
plt.plot(path_array[:,0],path_array[:,1],'b')
plt.show()
print('pause')

'''

# for step in np.arange(0,L2[0]*r+0.5,0.5):
#     e=nex_pos_i(dub,step,L2,'LSL',r)
#     path.append(e)

# path_array = np.array(path)
# plt.plot(p1[0],p1[1],'ro')
# plt.plot(p2[0],p2[1],'ro')
# plt.plot(path_array[:,0],path_array[:,1],'b')
# plt.show()
# print('pause')








