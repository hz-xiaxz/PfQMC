import numpy as np
from tenpy.networks.mps import MPS
from tenpy.models.hubbard import FermiHubbardModel
from tenpy.algorithms import dmrg

# 1. 设置模型参数
# TeNPy 的模型通常通过一个字典来配置，这非常方便。
# 这个字典的键名对应于你贴出文档中的 'Options'。
L = 12  # 系统长度
model_params = {
    'L': L,
    't': 1.0,
    'U': 4.0,
    'mu': 0.0,  # 化学势，通常用于控制粒子数，这里我们用量子数控制
    'cons_N': 'N',      # 开启粒子数守恒
    'cons_Sz': 'Sz',    # 开启总自旋 Z 分量守恒
    'bc_MPS': 'finite'  # 使用有限尺寸的 MPS 边界条件（对于 DMRG 是标准的）
    # 注意：'bc_x' 或 'bc_y' 可以用来设置模型的周期性边界条件，
    # 但 DMRG 在周期性边界上的收敛性较差。对于基态计算，开边界是首选。
}

print("TeNPy 程序：计算一维 Fermi-Hubbard 模型")
print("模型参数:")
for key, val in model_params.items():
    print(f"  {key}: {val}")
print("-" * 30)

# 2. 初始化模型
# TeNPy 会根据参数自动构建 MPO 形式的哈密顿量 H
# 它会自动处理 Jordan-Wigner 变换，我们无需关心。
M = FermiHubbardModel(model_params)

# 3. 创建初始波函数 (MPS)
# 我们需要一个处于正确量子数扇区的初始态。
# 这里我们选择半满（N=L=12 个电子）且总自旋为零 (Sz=0) 的扇区。
# 一个简单的初始态是交替的自旋构型 |↑,↓,↑,↓,...>
initial_state_list = ["up", "down"] * (L // 2)
# `from_product_state` 是一个强大的函数，可以从直积态创建 MPS。
# M.lat.mps_sites() 提供了格点的本地希尔伯特空间信息。
psi = MPS.from_product_state(M.lat.mps_sites(), initial_state_list, bc=M.lat.bc_MPS)


# 4. 设置并运行 DMRG 算法
# DMRG 的参数也通过一个字典来配置。
dmrg_params = {
    'mixer': True,  # 启用子空间扩展（mixer），帮助算法跳出局部最优解
    'max_sweeps': 10,
    'trunc_params': {
        'chi_max': 200,      # 每个 sweep 的最大键维度
        'svd_min': 1.e-10,   # 奇异值分解的截断阈值
    },
    'combine': True # 优化计算效率
}

print("开始进行 DMRG 计算...")
# 初始化 DMRG 引擎
eng = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
# 运行 DMRG
# `eng.run()` 返回基态能量和最终的基态波函数 psi
E, psi = eng.run()
print("DMRG 计算完成!")
print("-" * 30)

# 5. 打印和分析结果
print(f"计算得到的基态能量 E = {E:.10f}")

# (可选) 计算其他物理量
# 例如，计算每个格点的双占据数 <n_up * n_down>
# 'Nupdn' 是在 FermiHubbardModel 中预定义好的算符名
# 兼容旧版本的写法
# 使用 Python 的列表推导式来循环所有格点
# L 是我们定义的系统长度
double_occupancy = [psi.expectation_value('Nupdn', i) for i in range(L)]
avg_double_occupancy = np.mean(double_occupancy)

# 计算每个格点的 Sz
sz_values = psi.expect_onsite_g('Sz')
avg_sz = np.mean(sz_values)

print(f"平均双占据数 <Nup*Ndn> = {avg_double_occupancy:.6f}")
print(f"平均自旋 <Sz> = {avg_sz:.6f} (验证 Sz=0 的态)")