import numpy as np
import math
import itertools


from partition_function import computeRegionIdx,computeRegionCenters

def define_spec_region(allCenters, sets, partition, borderOutside=False):
    '''
    Return the indices of regions associated with the unique centers.

    Parameters
    ----------
    allCenters : List, e.g. 
    Ab.partition['R']['c_tuple'][(50.0, -2.98451, -2.98451, -140.0, -140.0, -2.98451)]=0
    center的坐标和在整体坐标中的序号
        List of the center coordinates for all regions.
    sets : List, 希望映射的区域
        List of sets to apply function for
    partition : Dict，之前已经得到的划分区域
        Partition dictionary.
    Returns
    -------
    list
        List of unique center points.

    '''
    
    delta = 1e-5
    
    if sets is None:
        return [], [], set()
    else:
    # set中每一个元素都描述了一个区域块
        points = [None] * len(sets)
        slices = {'min': [None] * len(sets), 'max': [None] * len(sets)}
        index_tuples = set()

        # Convert regions to all individual points (centers of regions)
        for i,set_boundary in enumerate(sets):#选取每一个单独的区域进行操作
        #若为str类型，那么说明该维度没有约束，因此用boundary即可, 如下表达每一个区域的空间部分，将'all'改成boundary的内容
            set_boundary = np.array([
                            S if type(S) != str else partition['boundary'][j] 
                            for j,S in enumerate(set_boundary) ])
            #向内缩小原始区域，来处理区域边界问题，只要足够小，那么就可以处理边界冲突的问题
            # Increase by small margin to avoid issues on region boundaries
            # print(set_boundary[:,[0]])
            set_boundary = np.hstack((set_boundary[:,[0]] + delta, 
                                      set_boundary[:,[1]] - delta))
            
            #区域的顶点坐标，这些顶点组成了一个凸多面体
            vertices = np.array(list(itertools.product(*set_boundary)))

            #根据给定的多面体顶点，给出包含着离散区域的序号，只保留留在序号范围内的值，
            _, indices_nonneg = computeRegionIdx(vertices, partition, borderOutside)
            #分别求最大和最小的序号，相当于矩形的两个顶点，i对应与维度序号
            slices['min'][i] = indices_nonneg.min(axis=0)
            slices['max'][i] = indices_nonneg.max(axis=0)

            # Create slices
            #分别利用每一个维度的最大和最小的两个顶点，生成其包围的所有区域的点，也就是原始区域的根据划分区域的包裹区域的center点
            # （注意，这些点都是顶点）
            index_tuples.update(set(itertools.product(*map(range, slices['min'][i], slices['max'][i]+1))))

            # Define iterator
            if borderOutside: # 内缩感兴趣区域
                M = map(np.arange, set_boundary[:,0]+delta, set_boundary[:,1]-delta, partition['width']/2)
            else:
                M = map(np.arange, set_boundary[:,0], set_boundary[:,1], partition['width']/2)

            # Retreive list of points and convert to array 重新组合原始边界的叉乘组合
            points[i] = np.array(list(itertools.product(*M)))

        points = np.concatenate(points)

        # Compute all centers of regions associated with points，计算上述边界点对应的center点，也就是顶点所在的区域的center
        centers = computeRegionCenters(points, partition)

        # Filter to only keep unique centers，按照行去除相同元素，如果同一个轴上的区域过小，那么经过抽象以后就会归结为同一个center
        # 此时，就会产生重复的points
        centers_unique = np.unique(centers, axis=0)

        #将所有坐标，映射成坐标的整体序号
        # for c in centers_unique:
        #     print(c)
        #     if tuple(c) in allCenters:
        #         print(c)
        states = [allCenters[tuple(c)] for c in centers_unique 
                           if tuple(c) in allCenters]
        
        # print('states:',len(states))
        # print('index_tuples:',len(index_tuples))
        if len(states) != len(index_tuples):
            print('ERROR: lengths of goal and goal_idx lists are not the same.')
            assert False

        # Return the ID's of regions associated with the unique centers 
        # state：绝对状态序号；slices：每个维度的上下值；index_tuples：边界center值        
        return states, slices, index_tuples
