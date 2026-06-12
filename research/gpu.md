
# 关于 $8 \times 4$ Hubbard Model 高精度模拟硬件配置的可行性分析报告


### 1. 摘要 (Executive Summary)

本报告旨在评估在 **Neural Quantum States (NQS)** 研究中，针对 **$8 \times 4$ Hubbard Model** 进行高精度基态模拟所需的显存门槛。

尽管我们计划采用先进的 **MinSR (Minimum-Norm Stochastic Reconfiguration)** 算法将显存复杂度从 $O(N_p^2)$ 降低至 $O(N_s^2)$，但基于 SOTA 文献（Phys. Rev. X / arXiv 2025）的标准配置，我们需要 **16,384** 的大采样数 ($N_s$) 来保证收敛精度。

**测算结论：** 在双精度复数 (Complex128) 环境下，16GB 显存的占用率将达到 **96% 以上**，处于极度危险的边缘，极易导致 OOM (内存溢出) 且无法进行网络扩展。**建议配置 24GB 显存 (RTX 5090)** 以确保科研项目的顺利进行。

---

### 2. 模拟参数设定 (Simulation Baseline)

为了复现顶刊级别的精度，我们设定以下不可妥协的物理参数：

*   **物理系统：** 2D Hubbard Model ($8 \times 4 = 32$ sites), 1/8 Doping.
*   **计算精度：** **Complex128** (双精度复数，16 Bytes/element)。
    *   *注：费米子符号问题要求波函数必须为复数，且化学精度要求双精度。*
*   **神经网络：** ResNet/GCNN 架构 (约 **150,000** 参数)。
*   **采样数 (Batch Size, $N_s$)：** **16,384**。
    *   *注：参考 SOTA 论文，只有达到此量级，VMC 的统计方差才能降至 $10^{-4}$ 量级，SR/MinSR 优化器才能稳定收敛。*

---

### 3. 显存开销详细核算 (Memory Breakdown)

即便使用 MinSR 算法避免了存储巨大的雅可比矩阵 ($N_s \times N_p$)，显存依然面临三大刚性支出的压力。

#### **A. 神经网络前向传播激活值 (Activations)**
无论使用何种优化算法，JAX 必须在显存中存储前向传播的中间结果以计算梯度 (Backpropagation)。
*   **估算：** 基于 16 层深、32 通道的 Deep GCNN。
*   **计算：** $N_s \times N_{sites} \times C \times L \times 16 \text{ Bytes}$
    $$ 16,384 \times 32 \times 32 \times 16 \times 16 \approx 4.3 \text{ GB} $$
*   **修正：** 考虑到 Attention Map (若涉及 Transformer) 及 JAX 的 Workspace，此项常驻显存约为 **6.0 GB**。

#### **B. MinSR 核心开销：Gram 矩阵 ($T T^\dagger$)**
MinSR 的核心是将优化问题转移到样本空间，需要构建并求解一个 $N_s \times N_s$ 的厄米矩阵 (Gram Matrix)。
*   **矩阵大小：** $N_s \times N_s$
*   **计算：** $16,384 \times 16,384 \times 16 \text{ Bytes}$
    $$ 268,435,456 \times 16 \text{ Bytes} \approx \mathbf{4.3 \text{ GB}} $$
*   *分析：* 这是 MinSR 算法带来的刚性开销，随着采样数 $N_s$ 的平方增长。

#### **C. 隐形开销：VJP 计算与系统占用**
*   **VJP Workspace：** MinSR 需要计算 $O_{loc} v$ (Jacobian-Vector Product)。虽然不需要存整个 Jacobian，但 JAX/XLA 在执行这一步大规模并行计算时，会预分配 **2.0 ~ 3.0 GB** 的临时显存。
*   **系统与框架：** 操作系统显示输出 + Python 解释器 + CUDA Kernels $\approx$ **2.0 GB**。

---

### 4. 风险评估：16GB vs 24GB

我们将上述开销汇总，对比两种硬件方案的可行性。

#### **方案一：RTX 5080 Laptop (16GB 显存)**

$$ \text{总需求} \approx \underbrace{2.0}_{\text{Sys}} + \underbrace{6.0}_{\text{Activations}} + \underbrace{4.3}_{\text{MinSR Matrix}} + \underbrace{3.0}_{\text{VJP/XLA}} = \mathbf{15.3 \text{ GB}} $$

*   **剩余空间：** **< 0.7 GB**。
*   **风险分析：**
    1.  **显存碎片 (Fragmentation)：** JAX 的内存分配器在 95% 以上占用率时，极易因为找不到连续内存块而报错 OOM。
    2.  **锁死扩展性：** 我们被死死限制在 150K 参数和 1.6万样本。如果为了更高精度想尝试 **PsiFormer (Attention)** 或增加网络宽度，显存直接溢出。
    3.  **调试困难：** 无法在 GPU 上同时运行调试工具或 Jupyter Notebook。

#### **方案二：RTX 5090 Laptop (24GB 显存)**

*   **剩余空间：** $24.0 - 15.3 = \mathbf{8.7 \text{ GB}}$。
*   **可行性分析：**
    1.  **安全冗余：** 拥有近 9GB 的动态空间，可以从容应对 JAX 的内存波动。
    2.  **科研潜力：** 允许我们将参数量翻倍至 **300K**，或引入 **Transformer** 架构（Attention 矩阵开销大），甚至将采样数 $N_s$ 进一步提升至 **20,000** 以追求极致精度。
    3.  **计算效率：** 充足的显存允许更大的 Chunk Size，减少 CPU-GPU 数据通信，训练速度显著提升。

---

### 5. 最终结论 (Conclusion)

虽然 **MinSR 算法** 成功规避了存储 36GB 雅可比矩阵的理论瓶颈，但在 **SOTA 级别的采样规模 ($N_s=16,384$)** 下，其自身的 Gram 矩阵开销配合神经网络的激活值，依然将 **16GB 显存** 逼入了死角。

*   **16GB (RTX 5080)：** 理论上刚好塞满，但在工程实践中极大概率无法稳定运行（OOM），且完全丧失了模型改进的空间。这是一个**“能跑代码，但出不了成果”**的配置。
*   **24GB (RTX 5090)：** 是能够稳定复现 $8 \times 4$ Hubbard Model 高精度结果的**最低物理门槛**。

**建议：** 鉴于目前有国家补贴政策，建议利用预算优先保障 **24GB 显存**，以确保研究生阶段科研任务的顺利完成。