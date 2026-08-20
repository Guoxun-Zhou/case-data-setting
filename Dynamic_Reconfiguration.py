# -*- coding: utf-8 -*-
"""

@author: zgx13
"""

from gurobipy import *
import numpy as np
from tupledictToarray import *
import math
'''
     （1）功能解释：
     本函数旨在实现考虑柔性负荷可调能力上的切负荷优化决策
     （2）输入dict
     dict_['branch_num'] #支路个数
     dict_['branch_innode'] #支路输入节点
     dict_['branch_outnode'] #支路流出节点
     dict_['branch_r'] #支路电阻
     dict_['branch_x'] #支路电抗
     dict_['node_num'] #节点个数
     dict_['node_P'] #节点有功
     dict_['node_Q'] #节点无功
     dict_['break_branch'] #故障支路
     dict_['break_branch_num'] #故障支路数
     dict_['node_DP'] #可调空间
     dict_['node_DP'] #电源点
'''

# %% 
class Load_trans_OPF():
    
    def __init__(self,dict_):
        self.branch_num = dict_['branch_num']
        self.branch_innode = dict_['branch_innode']
        self.branch_outnode = dict_['branch_outnode']
        self.branch_r = dict_['branch_r']
        self.branch_x = dict_['branch_x']
        self.node_num = dict_['node_num']
        self.node_P = dict_['node_P']
        self.node_Q = dict_['node_Q']
        self.node_PV_P = dict_['node_PV_P']
        self.T = dict_['T']
        # self.tttt = dict_['tttt']
        self.nosw = dict_['nosw']
        self.node_DP = dict_['node_DP']
        self.node_DQ = dict_['node_DQ']
        self.power_source = dict_['power_source']
        # self.balance_n = dict_['power_source']
        self.branch_Pmax = dict_['branch_Pmax']
        self.branch_Qmax = dict_['branch_Qmax']
        self.Uref = dict_['Umax']
        self.Us_UB = dict_['Umax']
        self.Us_LB = dict_['Umin']
        self.M = 100
# %%%% 构建模型求解
    def model_solve(self):
        
        # %%%%%% 创建模型
        model = Model('opf')
        # %%%%%% 创建开关变量
        z = model.addVars(self.branch_num,self.T,vtype = GRB.BINARY, name = 'z')
        # %%%%% 动作次数约束相关变量
        # wz = model.addVars(self.branch_num, self.T-1,vtype = GRB.BINARY, name='wz')
        # N = model.addVars(self.T-1,lb=0,vtype = GRB.CONTINUOUS, name = 'N')
        # w = model.addVars(self.T-1,vtype = GRB.BINARY, name = 'w')
        # S = model.addVars(self.T-1,vtype = GRB.BINARY, name = 'S')
        # %%%%%% 创建辅助的开关变量（用作分析辐射状）
        fz = model.addVars(self.branch_num,self.T,vtype = GRB.CONTINUOUS,lb = 0 ,ub = 1, name = 'fz')
        ff = model.addVars(self.branch_num,self.T,vtype = GRB.CONTINUOUS,lb = 0 ,ub = 1, name = 'ff')
        # %%%%%% 创建节点负荷变量
        P = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'P')
        Q = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'Q')
        # %%%%%% 创建节点切负荷变量
        Plr = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'Plr')
        Qlr = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'Qlr')
        # %%%%%% 创建光伏节点变量
        P_PV = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'P_PV')
        # %%%%%% 创建节点电压变量 Us = U^2
        Us = model.addVars(self.node_num,self.T,vtype = GRB.CONTINUOUS, lb = 0,name='Us')
        # %%%%%% 创建支路潮流变量
        BP = model.addVars(self.branch_num,self.T,vtype = GRB.CONTINUOUS,lb = -GRB.INFINITY, name = 'BP')
        BQ = model.addVars(self.branch_num,self.T,vtype = GRB.CONTINUOUS,lb = -GRB.INFINITY, name = 'BQ')
        # %%%%%% 创建支路电流变量 Is = I^2
        Is = model.addVars(self.branch_num,self.T,vtype = GRB.CONTINUOUS, lb = 0,name='Is')
        # %%% 整体开关动作约束
        # for t in range(self.T-1):
        #     model.addConstrs(wz[n,t] >= z[n,t+1]-z[n,t] for n in range(self.branch_num))
        #     model.addConstrs(wz[n,t] >= -z[n,t+1]+z[n,t] for n in range(self.branch_num))
        #     model.addConstr(N[t] == quicksum(wz[n,t] for n in range(self.branch_num)))
        #     model.addConstr(N[t] <= self.M*w[t])
        #     model.addConstr(N[t] >= 1 - self.M*(1-w[t]))
        #     model.addConstr(S[t] >= 1-self.M*(1-w[t]))
        #     model.addConstr(S[t] <= self.M*w[t])
        
        # model.addConstr(quicksum(S[t] for t in range(self.T-1)) <= 4)
        
        # %% 开关初始状态约束
        # model.addConstrs((z[n,t] == 1 for n in range(self.branch_num)  if n in range(self.node_num-3) for t in range(self.T)) ,name='开关初始状态')
        
        
        # %%% 切负荷约束 and 光伏出力约束
        for t in range(self.T):
            for j in range(self.node_num):
                # 光伏出力约束
                model.addConstr(P_PV[j,t] <= 0)
                model.addConstr(P_PV[j,t] >= self.node_PV_P[j,t])
                if j not in self.power_source:
                    # 切负荷约束
                    model.addConstr(Plr[j,t] >= 0)
                    model.addConstr(Plr[j,t] <= 0.1*self.node_P[j,t])
                    model.addConstr(Qlr[j,t] == Plr[j,t]*self.node_Q[j,t]/self.node_P[j,t])
                else:
                    model.addConstr(Plr[j,t] == 0)
                    model.addConstr(Qlr[j,t] == 0)
        #%% 节点负荷约束
        model.addConstrs((P[j,t] == self.node_P[j,t] - Plr[j,t] + P_PV[j,t] for j in range(self.node_num) for t in range(self.T)  if j not in self.power_source),name='节点有功功率约束'  )
        model.addConstrs((Q[j,t] == self.node_Q[j,t] - Qlr[j,t] for j in range(self.node_num) for t in range(self.T)  if j not in self.power_source),name='节点无功功率约束'  )
        # %% 电流约束
        model.addConstrs(Is[i,t]<=10*z[i,t] for i in range(self.branch_num) for t in range(self.T))
        # %% 分时段开关限制
        # a = [4,4,4,4,2,2,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,1,1]  #FCM
        a = [1,1,1,1,1,1,0,0,0,0,3,3,3,3,3,2,2,2,2,4,4,4,4,4]  # 决策
        classification_0,classification_1,classification_2,classification_3,classification_4,classification_5,classification_6 = [],[],[],[],[],[],[]
        for i in range(24):
            if a[i] == 0:
                classification_0.append(i)
            elif a[i] == 1:
                classification_1.append(i)
            elif a[i] == 2:
                classification_2.append(i)
            elif a[i] == 3:
                classification_3.append(i)
            elif a[i] == 4:
                classification_4.append(i)
            elif a[i] == 5:
                classification_5.append(i)
            elif a[i] == 6:
                classification_6.append(i)
            else:
                pass
        for n in range(self.branch_num):
            model.addConstrs(z[n,classification_0[0]] == z[n,t] for t in classification_0)
            model.addConstrs(z[n,classification_1[0]] == z[n,t] for t in classification_1)
            model.addConstrs(z[n,classification_2[0]] == z[n,t] for t in classification_2)
            model.addConstrs(z[n,classification_3[0]] == z[n,t] for t in classification_3)
            model.addConstrs(z[n,classification_4[0]] == z[n,t] for t in classification_4)
            model.addConstrs(z[n,classification_5[0]] == z[n,t] for t in classification_5)
            model.addConstrs(z[n,classification_6[0]] == z[n,t] for t in classification_6)
        # %% 常闭开关线路约束
        for t in range(self.T):
            model.addConstrs(z[n,t] == 1 for n in self.nosw[t] )
        # %% 功率上下界约束
        for t in range(self.T):
            for n in range(self.node_num):
                if n in self.power_source:
                    model.addConstr((self.node_P[n,t] - Plr[n,t] + P_PV[n,t] - self.node_DP[n] - P[n,t] <= 0.0  ),name = '节点有功下界约束')
                    model.addConstr((self.node_Q[n,t] - Qlr[n,t] - self.node_DQ[n] - Q[n,t] <= 0.0  ),name = '节点无功下界约束')
                    model.addConstr((self.node_P[n,t] - Plr[n,t] + P_PV[n,t] - P[n,t] >= 0  ),name = '节点有功上界约束')
                    model.addConstr((self.node_Q[n,t] - Qlr[n,t] - Q[n,t] + 10000>= 0  ),name = '节点无功上界约束')
                else:
                    model.addConstr((self.node_P[n,t] - Plr[n,t] + P_PV[n,t] - self.node_DP[n] - P[n,t] <= 0.0  ),name = '节点有功下界约束')
                    model.addConstr((self.node_Q[n,t] - Qlr[n,t] - self.node_DQ[n] - Q[n,t] <= 0.0  ),name = '节点无功下界约束')
                    model.addConstr((self.node_P[n,t] - Plr[n,t] + P_PV[n,t] - P[n,t] >= 0  ),name = '节点有功上界约束')
                    model.addConstr((self.node_Q[n,t] - Qlr[n,t] - Q[n,t] >= 0  ),name = '节点无功上界约束')
        # %%%%%% 创建辐射状连通性约束（及部分有功无功约束）
            # model.addConstr(quicksum(z[n,t] for n in range(self.branch_num)) == self.node_num - 4)
            model.addConstrs(fz[b,t] + ff[b,t] == z[b,t]  for b in range(self.branch_num)) 
            mid_ffs = 0.0
            mid_fzs = 0.0
            # node_branchin = [] # 注入节点的支路编号 
            # node_branchout = [] # 节点流出至其他支路编号
            for n in range(self.node_num):
                mid_branch_in = []
                mid_branch_out = []
                if n not in self.power_source:
                    mid_ff = 0.0
                    mid_fz = 0.0
                    for b in range(self.branch_num):
                        if self.branch_innode[b] == n:
                            mid_ff += ff[b,t]
                            mid_branch_out.append(b)
                        if self.branch_outnode[b] == n:
                            mid_fz += fz[b,t]
                            mid_branch_in.append(b)
                    model.addConstr(mid_fz + mid_ff <= (self.node_num - len(self.power_source)) / (self.node_num - len(self.power_source) + 1)  , name='3b') # 辐射状约束（非上层电网注入节点）
                else:
                    for b in range(self.branch_num):
                        if self.branch_innode[b] == n:
                            mid_ffs += ff[b,t]
                            mid_branch_out.append(b)
                        if self.branch_outnode[b] == n:
                            mid_fzs += fz[b,t]
                            mid_branch_in.append(b) 
                model.addConstr((quicksum(BP[b,t] for b in mid_branch_out) + P[n,t] == quicksum(BP[b,t] - self.branch_r[b]*Is[b,t] for b in mid_branch_in)),name='有功平衡约束') #有功平衡约束
                model.addConstr((quicksum(BQ[b,t] for b in mid_branch_out) + Q[n,t] == quicksum(BQ[b,t] - self.branch_x[b]*Is[b,t] for b in mid_branch_in)),name='无功平衡约束') #无功平衡约束
            model.addConstr(mid_fzs + mid_ffs <= (self.node_num - len(self.power_source)) / (self.node_num - len(self.power_source)+ 1)  , name='3c')# 辐射状约束（上层电网注入节点）
        
        # %%%%%%%% 节点电压约束
            model.addConstrs( 0.0 + self.M*(1-z[i,t]) >= -Us[self.branch_outnode[i],t] + Us[self.branch_innode[i],t] - 2*(self.branch_r[i]*BP[i,t]+self.branch_x[i]*BQ[i,t])
                             + (self.branch_r[i]**2 + self.branch_x[i]**2)*Is[i,t]  for i in range(self.branch_num))
            model.addConstrs( 0.0 - self.M*(1-z[i,t]) <= -Us[self.branch_outnode[i],t] + Us[self.branch_innode[i],t] - 2*(self.branch_r[i]*BP[i,t]+self.branch_x[i]*BQ[i,t])
                             + (self.branch_r[i]**2 + self.branch_x[i]**2)*Is[i,t]  for i in range(self.branch_num))
            model.addConstrs(Us[i,t] == self.Uref for i in self.power_source)
            
            model.addConstrs(Us[i,t] <= self.Us_UB for i in range(self.node_num))
            model.addConstrs(Us[i,t] >= self.Us_LB for i in range(self.node_num))
        
        # %%%%%% 额外约束
            model.addConstrs(4*BP[i,t]*BP[i,t] + 4*BQ[i,t]*BQ[i,t] + 
                              (Is[i,t]-Us[self.branch_innode[i],t])*(Is[i,t]-Us[self.branch_innode[i],t]) 
                              <= (Is[i,t] + Us[self.branch_innode[i],t]) * (Is[i,t]+Us[self.branch_innode[i],t])
                                  for i in range(self.branch_num))

        # %%%%%% 模型目标
        model.setObjective(quicksum(1000*Is[n,t]*self.branch_r[n] for n in range(self.branch_num) for t in range(self.T))
                          +7*1000*quicksum(quicksum((Is[k,t] -quicksum(Is[i,t] for i in range(self.branch_num))/int(self.branch_num))**2  for k in range(self.branch_num))/ int(self.branch_num) for t in range(self.T))
                            +quicksum(P_PV[j,t]-self.node_PV_P[j,t] for j in range(self.node_num) for t in range(self.T))
                            # +0.000001*quicksum(wz[n,t] for n in range(self.branch_num) for t in range(self.T-1)) 
                           ,GRB.MINIMIZE)
        model.Params.NonConvex = 2
        model.write('cut.lp')
        model.optimize()
        # model.computeIIS()
        # model.write('model.ilp')
        result_ = {}
        result_['obj_'] = model.getObjective().getValue()
        # obj_ = model.getObjective().getValue()
        
        result_['z'] = double_var(z,self.branch_num,self.T)
        result_['P'] = double_var(P,self.node_num,self.T)
        result_['Q'] = double_var(Q,self.node_num,self.T)
        result_['Plr'] = double_var(Plr,self.node_num,self.T)
        result_['P_PV'] = double_var(P_PV,self.node_num,self.T)
        result_['Qlr'] = double_var(Qlr,self.node_num,self.T)
        result_['P_branch'] = double_var(BP,self.branch_num,self.T)
        result_['Q_branch'] = double_var(BQ,self.branch_num,self.T)
        result_['Us'] = double_var(Us,self.node_num,self.T)
        result_['Is'] = double_var(Is,self.branch_num,self.T)
        result_['Is_s'] = np.sqrt(result_['Is'])

        result_['ff'] = double_var(ff,self.branch_num,self.T)
        result_['fz'] = double_var(fz,self.branch_num,self.T)
        # result_['x2'] = double_var(self.x2_all,self.n2,1)
        # obj_L.append(result_['obj_'])
        print('----------------成功运行------------------')
        # print('目标网损成本为：',quicksum(1000*result_['Is'][n,t] for n in range(self.branch_num) for t in range(self.T)))
        print('目标负载均衡指数为：',7000*quicksum(quicksum((result_['Is'][k,t] -quicksum(result_['Is'][i,t] for i in range(self.branch_num))/int(self.branch_num))**2  for k in range(self.branch_num))/ int(self.branch_num) for t in range(self.T)))
        print('网损成本为：',quicksum(1000*result_['Is'][n,t]*self.branch_r[n] for n in range(self.branch_num) for t in range(self.T)))
        # print('负载均衡指数为：',7000*quicksum(quicksum((result_['Is_s'][k,t] -quicksum(result_['Is_s'][i,t] for i in range(self.branch_num))/int(self.branch_num))**2  for k in range(self.branch_num))/ int(self.branch_num) for t in range(self.T)))
        loss_all = []
        for t in range(self.T):
            loss_all_flag = sum(1000*result_['Is'][n,t]*self.branch_r[n] for n in range(self.branch_num))
            loss_all.append(loss_all_flag)
            loss_all_flag = 0
        print('各时刻网损结果为：',loss_all)
        
        load_piancha = []
        for t in range(self.T):
            load_list = [-1000*result_['P'][0,t],-1000*result_['P'][1,t],-1000*result_['P'][2,t]]
            load_flag = np.std(load_list)
            load_piancha.append(load_flag)
            load_flag = 0
            load_list = []
        print('各时刻负载均衡结果为：',load_piancha)
        return result_
        
        
