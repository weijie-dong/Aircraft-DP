import numpy as np
import math
import itertools        


# 可获得每个region的center和整体的序号
def define_partition(dim, nrPerDim, regionWidth, origin):
    '''
    Define the partitions object `partitions` based on given settings.

    Parameters
    ----------
    dim : int
        Dimension of the state (of the LTI system).
    nrPerDim : list
        每一维度需要划分的区域个数
        List of integers, where each value is the number of regions in that 
        dimension.
    regionWidth : list
        每一维度区域的宽度
        Width of the regions in every dimension.
    origin : list
        每一维度原点位置
        Coordinates of the origin of the continuous state space. 

    Returns
    -------
    partition : dict
        Dictionary containing the info regarding partisions.

    '''
    print('act partition process')
    regionWidth = np.array(regionWidth)
    origin      = np.array(origin)
    elemVector   = dict()
    for i in range(dim):#遍历每一个维度
        
        elemVector[i] = np.linspace(-(nrPerDim[i]-1)/2, 
                                     (nrPerDim[i]-1)/2,
                                     int(nrPerDim[i]))
        #每个维度的划分的每个格子的中点，以整体区域的中点作为原点，而非原始系统的原点      
    widthArrays = [[x*regionWidth[i] for x in elemVector[i]] 
                                              for i in range(dim)]
    #为每个格子的中点进行标号，到序号的映射
    idxArrays = [range(len(arr)) for arr in widthArrays]
     #划分的region的总数，prod求乘积，获得region的数目
    nr_regions = np.prod(nrPerDim)
    #区域保存的初始化
    partition = {'center': np.zeros((nr_regions, dim), dtype=float), 
                 'idx': {},
                 'idx_inv': {},
                 'c_tuple': {}}
    #包含center左下和右上两个对角的点
    partition['low'] = np.zeros((nr_regions, dim), dtype=float)
    partition['upp'] = np.zeros((nr_regions, dim), dtype=float)

    count=0
    # 遍历格子（中点，对应序号）
    # print('widthArrays:',widthArrays,'\n')
    print('idxArrays:',idxArrays,'\n')
    for i,(pos,idx) in enumerate(zip(itertools.product(*widthArrays),
                                     itertools.product(*idxArrays))):
        # if i//200==0:
        #     print('compute each dim')
        # count+=1
        # if count>20: break
        center = np.array(pos) + origin
        #保留小数位数
        dec = 5
        partition['center'][i] = np.round(center, decimals=dec)
        partition['c_tuple'][tuple(np.round(center, decimals=dec))] = i #每个cernter坐标到整体序号的映射
        partition['low'][i] = np.round(center - regionWidth/2, #m每个center所在的上下界
                                        decimals=dec)
        partition['upp'][i] = np.round(center + regionWidth/2, 
                                        decimals=dec)
        partition['idx'][tuple(idx)] = i #每个维度上的坐标到整体序号的映射
        partition['idx_inv'][i] = tuple(idx)#整体序号到每个维度上的坐标的映射

    return partition


def computeRegionIdx(points, partition, borderOutside=False):
    '''
    Function to compute the indices of the regions that a list of points belong

    Parameters
    ----------
    points : 2D Numpy array，记录了一个多面体的所有顶点
        Array, with every row being a point to determine the center point for.
    partition : dict
        Dictionary of the partition.

    Returns
    -------
    2D Numpy array
        Array, with every row being the indices.这个indices是每一维上隶属的状态的序号

    '''
    
    # Shift the points to account for a non-zero origin
    #此处先求距离origin的距离，然后整体向右平移partition['width']*partition['number']/2距离
    #相当于将最左侧的边界点，挪到原点处后的新的坐标，这是一个绝对坐标，为了计算后续属于该维度的第几个分区
    pointsZero = points - partition['origin'] + \
                    partition['width']*partition['number']/2
    #平移过后再除以partition['width']，即可获得属于那个维度的第几个区域，这对应于每个维度的第几个划分
    indices = (pointsZero // partition['width']).astype(int)
    # Reduce index by one if it is exactly on the border
    #borderOutside=False 则即使落在边界，也不去除indices
    indices -= ((pointsZero % partition['width'] == 0).T * borderOutside).T
    #超出边界的点，就直接保留边界，有可能存在负值，但好像不会超出最大边界，因为是向内收缩
    indices_nonneg = np.minimum(np.maximum(0, indices), 
                                np.array(partition['number'])-1).astype(int)
    
    return indices, indices_nonneg

def computeRegionCenters(points, partition):
    '''
    Function to compute to which region (center) a list of points belong

    Parameters
    ----------
    points : 2D Numpy array
        Array, with every row being a point to determine the center point for.
    partition : dict
        Dictionary of the partition.

    Returns
    -------
    2D Numpy array
        Array, with every row being the center coordinate of that row of the
        input array.

    '''
    
    # Check if 'points' is a vector or matrix，将一维的全部转成二维的结构
    if len(np.shape(points)) == 1:
        points = np.reshape(points, (1,len(points)))

    # Retreive partition parameters，提取相关参数
    region_width = np.array(partition['width'])
    region_nrPerDim = partition['number']
    dim = len(region_width)

    # Boolean list per dimension if it has a region with the origin as center，origin是center中的一个的序列，
    # origin是center那么是False
    originCentered = [True if nr % 2 != 0 else False for nr in region_nrPerDim]

    # Initialize centers array
    centers = np.zeros(np.shape(points)) 

    # Shift the points to account for a non-zero origin
    originShift = np.array(partition['origin'] )
    pointsShift = points - originShift #计算points和中心点的距离

    for q in range(dim):
        # Compute the center coordinates of every shifted point
        if originCentered[q]:
            #若center不是origin，也就是说origin是边界
            centers[:,q] = ((pointsShift[:,q]+0.5*region_width[q]) //
                             region_width[q]) * region_width[q]
        
        else:
            # origin是center，因此先把坐标向左移半格
            centers[:,q] = (pointsShift[:,q] // region_width[q]) * \
                            region_width[q] + 0.5*region_width[q]
    # Add the origin again to obtain the absolute center coordinates，获得绝对坐标
    return np.round(centers + originShift, decimals=5)