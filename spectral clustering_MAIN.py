# -*- coding: utf-8 -*-
"""
@author: zgx13
"""
import time
import copy
import pandas as pd
import numpy as np
import Dynamic_Reconfiguration
import Dynamic_Reconfiguration3
import random
import Plot_Dynamic
import xlwt


def TO_input2(path,t,z):
    
    '''
    读取标准算例数据
    '''
    dict_ = {}
    # %% 标幺化参数设置
    V_base = 12.66 
    S_base = 5.68  # MVA 33:5.68/ 94:5.68
    R_base = V_base**2 / S_base
    # %% 构建联络线路的拓扑信息
    # %%%% 加载文件
    branch = pd.read_excel(path+'branch.xlsx',header=0,index_col=0).values
    # %%%% 支路数目
    dict_['branch_num'] = len(branch)
    # %%%% 注入流出节点
    dict_['branch_innode'] = (branch[:,0]).astype(int)
    dict_['branch_outnode'] = (branch[:,1]).astype(int)
    # %%%% 电阻电抗
    dict_['branch_r'] = branch[:,2] / R_base
    dict_['branch_x'] = branch[:,3] / R_base
    # %%%% 支路容量极限
    dict_['branch_Pmax'] = 4.3*np.ones(dict_['branch_num']) #后面需要改动
    dict_['branch_Qmax'] = 4.3*np.ones(dict_['branch_num']) # 后面需要改动
    
    # %%节点功率数据
    load = pd.read_excel(path + 'load.xlsx', header=0, index_col=0).values
    dict_['node_num'] = len(load)
    # dict_['node_P'] = load[:, 0] / S_base / 1000
    # dict_['node_Q'] = load[:, 1] / S_base / 1000
    load_P = pd.read_excel(path +'load_P.xlsx').values
    load_Q = pd.read_excel(path +'load_Q.xlsx').values
    dict_['node_P'] = load_P/S_base / 1000
    dict_['node_Q'] = load_Q/S_base / 1000
    # %% 光伏数据
    PV_P = pd.read_excel(path +'PV.xlsx').values
    dict_['node_PV_P'] = PV_P/S_base / 1000
    
    dict_['T'] = t
    # %%%% 电源节点集合
    dict_['power_source'] = np.arange(3).astype(int)# 33:[0] 94:0-10
    # %%%%PV节点集合
    dict_['PV_node'] = [9]
    # %%%% 调整联络线的容量约束
    dict_['branch_Pmax'][0] = 10 # 后面需要改动
    dict_['branch_Qmax'][0] = 10 # 后面需要改动
    
    dict_['Umax'] = 1.05**2
    dict_['Umin'] = 0.95**2
    # %%%%  节点可下调空间(可调整)
    dict_['node_DP'] = np.zeros(dict_['node_num'])
    dict_['node_DQ'] = np.zeros(dict_['node_num'])
    for i in dict_['power_source']:
        dict_['node_DP'][i] = 200 # 等效于平衡节点的输送容量大小为20
        dict_['node_DQ'][i] = 200 # 等效于平衡节点的输送容量大小为20
    
    # %% 设置开关常闭线路
    # Nosw = pd.read_excel(path +'stay_close.xlsx', header=None, index_col=None).values
    dict_['nosw'] = []
    dict_['nosw_num'] = len(dict_['nosw'])
    
    # %% 导入确定的结果
    dict_['zzz'] = z
    return dict_

    # %%%%
if __name__ == '__main__':
    time_start = time.time()
    # %% 求取矩阵
    path = 'trans_case/ieee16_1/' # 
    result_matrix = np.zeros([24,24])
    
    result_dict = {}
    # 预计算所有时刻的OS结果
    for i in range(24):
        z_z = []
        dict_2 = TO_input2(path, i, z_z)
        OS = Dynamic_Reconfiguration.Load_trans_OPF(dict_2)
        result2 = OS.model_solve()
        
        # print(result2['z'])
        # print(type(result2['z']))
        
        # 提取单时刻重构的解，保存在字典中
        rounded_array = result2['z'].round().astype(int)
        result_dict[i] = rounded_array
    
    result_matrix = np.zeros((24, 24))
    
    # 使用预计算的结果进行后续计算
    for i in range(24):
        z_z = result_dict[i]  # 直接使用预计算的结果
        
        for j in range(24):
            dict_3 = TO_input2(path, j, z_z)
            SS = Dynamic_Reconfiguration3.Load_trans_OPF(dict_3)
            result3 = SS.model_solve()
            if i > j:
                result_matrix[i,j] = (1+0.2*(i-j))*result3['obj_'] + 0.005*np.sum((result_dict[i] - result_dict[j]) ** 2)    #(1+0.02*(i-j))*
            elif i < j:
                result_matrix[i,j] = (1+0.2*(j-i))*result3['obj_'] + 0.005*np.sum((result_dict[i] - result_dict[j]) ** 2)    #(1+0.02*(j-i))*
            else:
                result_matrix[i,j] = result3['obj_']
                
            
            '''
            注：0.12*（i-j）和0.12*（j-i）是为了考虑时序性,取值是暂定的，可以根据实际进行修改
            '''
    result_matrix_1 = result_matrix.copy()
    # print('result_matrix_1的结果')
    # print(result_matrix_1)
    for i in range(24):
        for j in range(24):
            result_matrix[i][j] = 100*result_matrix_1[i][j]/result_matrix_1[j][j] - 99
    
    # print('--------最终的结果矩阵如下----------------')
    # print(result_matrix)
    
    
    # %% 对矩阵进行排序，计算相似度，得到分段
    step1 = result_matrix
    step2 = np.zeros([24,24])
    for i in range(24):
        for j in range(24):
            step2[i][j] = step1[i][j] - step1[j][j]
            
    # print(step2)
    print('-------------------------')
    step22 = np.zeros([24,24])
    flag = 0
    for i in range(24):
        for j in range(24):
            flag = step2[i][j] +step2[j][i]
            step22[i][j] = flag
    
    # print('step22:')
    # print(step22)
    #将相似度矩阵输出到excel中
    f = xlwt.Workbook('encoding = utf-8') #设置工作簿编码
    sheet1 = f.add_sheet('sheet1',cell_overwrite_ok=True) #创建sheet工作表
    # list1 = a #要写入的列表的值
    for j in range(24):
        for i in range(24):
            list1 = list(step22) #要写入的列表的值
            sheet1.write(j,i,list1[j][i]) #写入数据参数对应 行, 列, 值
    f.save('similarity_matrix.xls')#保存.xls到当前工作目录

    step3 = np.zeros([24,24])
    flag = 0
    for i in range(24):
        for j in range(24):
            flag = step2[i][j] +step2[j][i]
            if flag == 0:
                step3[i][j] = 50
            else:
                step3[i][j] = 100/flag
    # print('这里是STEP3')    
    # print(step3)
    


    from sklearn.cluster import SpectralClustering

    # 假设相似度矩阵为similarity_matrix
    # similarity_matrix是一个二维数组，表示样本之间的相似度

    # 创建谱聚类对象
    spectral_clustering = SpectralClustering(n_clusters=5, affinity='precomputed')

    # 进行聚类
    similarity_matrix = step3
    labels = spectral_clustering.fit_predict(similarity_matrix)

    # 输出每个样本的聚类标签
    print('分段结果为：',labels)
    
    # %% 将分段代入到
    
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('运行的时间为：',time_sum)
    # # 求解
    # LT = Dynamic_Reconfiguration2.Load_trans_OPF(dict_)
    # result = LT.model_solve()
    
    # 输出最后保存的开关结果
    print('最后保存的开关结果：',result_dict)
    
    # %%%% 画图
    
    # dict_plot = {}
    # dict_plot['path'] = path
    # dict_plot['node_num'] = dict_['node_num']
    # dict_plot['branch_num'] = dict_['branch_num']
    # dict_plot['branch_innode'] = dict_['branch_innode']
    # dict_plot['branch_outnode'] = dict_['branch_outnode']
    # dict_plot['T'] = dict_['T']
    # dict_plot['w'] = result['z']
    # dict_plot['G_size'] = (15,10)
    # Plot_Dynamic.plot_net(dict_plot)