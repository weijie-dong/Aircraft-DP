import numpy as np
import math
import itertools

from spec_region import define_spec_region
from partition_function import define_partition
from model_spec import AIRCRAFT, UAV_spec, Input_spec

class Abstraction(object):

    def __init__(self):
        self.model = AIRCRAFT()
        self.spec=UAV_spec()

    def partion_state(self):
        self.partition=dict()
        self.partition['state_variables']=self.model.state_variables
        # Determine origin of partition
        # number of region in each dim
        self.spec.partition['number'] = np.array(self.spec.partition['number'])
        # width of each region in each dim
        self.spec.partition['width'] = np.array([-1,1]) @ \
            self.spec.partition['boundary'].T / self.spec.partition['number']
        
        #find the origin point
        self.spec.partition['origin'] = 0.5 * np.ones(2) @ \
            self.spec.partition['boundary'].T
        #判断origin在每一维度是否是当前维度的中心，True表示origin不是center之一
        self.spec.partition['origin_centered'] = \
            [True if nr % 2 != 0 else False for nr in self.spec.partition['number']]
        
        print(' -- Width of regions in each dimension:',self.spec.partition['width'])

        # 可获得每个region的center和整体的序号
        self.partition['R'] = define_partition(self.model.n,
                            self.spec.partition['number'],
                            self.spec.partition['width'],
                            self.spec.partition['origin'])
        
        #获得region的数量
        self.partition['nr_regions'] = len(self.partition['R']['center'])

        print(' -- Number of regions:',self.partition['nr_regions'])

        # Create all combinations of n bits, to reflect combinations of all 
        # lower/upper bounds of partitions 2的n次方个值，所有0-1的组合
        bitCombinations = list(itertools.product([0, 1], 
                               repeat=self.model.n))
        # print(bitCombinations)
        bitRelation     = ['low','upp']
        # for i in range(self.partition['nr_regions']):
        #     for bitList in bitCombinations:
        #         for bitIndex,bit in enumerate(bitList):
        #             print(bitRelation[bit])
        #             print(self.partition['R'][bitRelation[bit]][i][bitIndex])

        # Calculate all corner points of every partition. Every partition has 
        # an upper and lower bounnd in every state (dimension). Hence, the 
        # every partition has 2^n corners, with n the number of states.
        # nr_regions是center的数量
        #bitIndex 表示第几个维度的值
        #bit表示up还是low
        # 计算的是每个center的周边的corners
        # 如何避免重复？
        self.partition['allCorners'] = np.array( [[[
            self.partition['R'][bitRelation[bit]][i][bitIndex] 
                for bitIndex,bit in enumerate(bitList)
            ] for bitList in bitCombinations 
            ] for i in range(self.partition['nr_regions'])
            ] )


        # 将目标区域进行包裹。goal记录其绝对状态序号
        # Determine goal regions
        self.partition['goal'], self.partition['goal_slices'], \
            self.partition['goal_idx'] = define_spec_region(
                allCenters = self.partition['R']['c_tuple'], 
                sets = self.spec.goal['ownship_t'],
                partition = self.spec.partition,
                borderOutside = True)
        
        print(' -- Number of goal regions:',len(self.partition['goal']))

        # 将危险区域进行包裹
        # 'critical'：绝对状态序号；'critical_slices'：每个维度的上下值；'critical_idx'：边界center值        
        self.partition['critical'], self.partition['critical_slices'], \
            self.partition['critical_idx'] = define_spec_region(
                allCenters = self.partition['R']['c_tuple'], 
                sets = self.spec.critical,
                partition = self.spec.partition,
                borderOutside = True)
        
        print(' -- Number of critical regions:',len(self.partition['critical']))

    def initialize_results(self):
        
        # Initialize results dictionaries
        self.results = dict()
        v = 30

        self.results['optimal_policy'] = np.zeros((v, self.partition['nr_regions']), dtype=int)


class Abstraction_u(object):

    def __init__(self):
        self.spec=Input_spec()

    def partion_state(self):
        self.partition=dict()
        # Determine origin of partition
        # number of region in each dim
        self.spec.partition['number'] = np.array(self.spec.partition['number'])
        # width of each region in each dim
        self.spec.partition['width'] = np.array([-1,1]) @ \
            self.spec.partition['boundary'].T / self.spec.partition['number']
        
        #find the origin point
        self.spec.partition['origin'] = 0.5 * np.ones(2) @ \
            self.spec.partition['boundary'].T
        #判断origin在每一维度是否是当前维度的中心，True表示origin不是center之一
        self.spec.partition['origin_centered'] = \
            [True if nr % 2 != 0 else False for nr in self.spec.partition['number']]
        
        print(' -- Width of regions in each dimension:',self.spec.partition['width'])
        dim=self.spec.partition['dim']
        # 可获得每个region的center和整体的序号
        self.partition['R'] = define_partition(dim,
                            self.spec.partition['number'],
                            self.spec.partition['width'],
                            self.spec.partition['origin'])
        
        #获得region的数量
        self.partition['nr_regions'] = len(self.partition['R']['center'])

        print(' -- Number of regions:',self.partition['nr_regions'])

        # Create all combinations of n bits, to reflect combinations of all 
        # lower/upper bounds of partitions 2的n次方个值，所有0-1的组合
        bitCombinations = list(itertools.product([0, 1], 
                               repeat=dim))
        # print(bitCombinations)
        bitRelation     = ['low','upp']
        # for i in range(self.partition['nr_regions']):
        #     for bitList in bitCombinations:
        #         for bitIndex,bit in enumerate(bitList):
        #             print(bitRelation[bit])
        #             print(self.partition['R'][bitRelation[bit]][i][bitIndex])

        # Calculate all corner points of every partition. Every partition has 
        # an upper and lower bounnd in every state (dimension). Hence, the 
        # every partition has 2^n corners, with n the number of states.
        # nr_regions是center的数量
        #bitIndex 表示第几个维度的值
        #bit表示up还是low
        # 计算的是每个center的周边的corners
        # 如何避免重复？
        self.partition['allCorners'] = np.array( [[[
            self.partition['R'][bitRelation[bit]][i][bitIndex] 
                for bitIndex,bit in enumerate(bitList)
            ] for bitList in bitCombinations 
            ] for i in range(self.partition['nr_regions'])
            ] )
        
    def initialize_results(self):
        
        # Initialize results dictionaries
        self.results = dict()
        v = 30

        self.results['optimal_policy'] = np.zeros((v, self.partition['nr_regions']), dtype=int)