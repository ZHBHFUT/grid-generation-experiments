# HyperGrid：二维双曲型 O 型结构网格生成器

独立 C++17 命令行程序。基于 Ricky Reusser 的 [grid-generation-experiments](https://github.com/rreusser/grid-generation-experiments) 中的双曲推进方程、周期差分及网格推进算法重新实现（现默认采用线性化隐式方法），输入由 NACA 参数改为任意简单闭合的逆时针二维多段线。原始算法说明：[Hyperbolic Grid Generation](https://rreusser.github.io/2016/02/07/hyperbolic-grid-generation/)。上游 `package.json` 声明 ISC 许可证；上游仓库未包含单独的 LICENSE 文件。本项目未复制原网页、npm 依赖或第三方库。

## 1. 编译

要求支持 C++17 的编译器与 CMake ≥ 3.16，无需第三方 C++ 库。Python 3 仅用于可选的回归测试。

Windows（Visual Studio Developer PowerShell 或已安装 CMake/Visual Studio 的 PowerShell）：

```powershell
cmake -S . -B build
cmake --build build --config Release
.\build\Release\hypergrid.exe --help
```

以上为 Visual Studio 多配置生成器。如果使用 Ninja 等单配置生成器：

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
.\build\hypergrid.exe --help
```

Linux：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/hypergrid --help
```

也可以直接编译单个源文件，例如 `g++ -std=c++17 -O3 src/main.cpp -o hypergrid`。

## 2. 使用方法

Windows PowerShell（Visual Studio 生成器）：

```powershell
.\build\Release\hypergrid.exe `
    --input examples/naca0012.csv `
    --output mesh.xyz `
    --first-height 1e-4 `
    --growth 1.15 `
    --layers 40
```

`--first-height` 与输入坐标使用同一长度单位。`--layers` 为**单元层数**，不含壁面这条节点线；输出共 `layers+1` 条节点线。

### 输入 CSV

```csv
x,y
1.0,0.0
0.8,0.04
0.4,0.06
0.0,0.0
0.4,-0.06
0.8,-0.04
1.0,0.0
```

- 仅接受 `x,y` 两列；允许大小写 `X,Y` 表头、空行、`#` 注释以及 UTF-8 BOM。
- 至少 5 个独立节点，按逆时针方向排列，简单闭合，不允许自交、回折重叠或相邻重复节点。
- 最后一行可以重复第一行，也可以省略；程序自动将轮廓闭合。
- 保留壁面节点，不自动重采样或改变其分布；对于极尖的角点、深凹轮廓、相邻点距骤变、极不均匀的点距，可能需要预处理。

### 参数

| 参数 | 默认 | 意义 |
|---|---|---|
| `--input PATH` | 必填 | CSV 文件 |
| `--output PATH` | `grid.xyz` | 输出 Plot3D ASCII 文件 |
| `--first-height H` | 必填 | 第一层物理高度（height 默认设置）或第一推进步长（original 默认设置） |
| `--growth Q` | 1.15 | 几何增长比，必须大于 0；1 表示等步长 |
| `--layers N` | 必填 | 单元层数 |
| `--mode height\|original` | height | 物理层高模式或原始 PDE 模式 |
| `--alpha A` | height: 0; original: 1 | 非负速度指数 |
| `--diffusion E` | 0.002 | 无量纲光顺系数，0 关闭光顺 |
| `--smoothing tangent\|full` | height: tangent; original: full | 仅保留切向光顺或完整二阶导数光顺 |
| `--smoothing-metric wall\|current` | wall | 切向光顺度量：初始壁面弧长差分，或当前层两侧实际点距的倒数差；current 仅支持 tangent 光顺 |
| `--smoothing-start H` | height: 0.005L; original: 0 | 名义推进距离小于该值时关闭光顺（输入长度单位） |
| `--smoothing-end H` | height: 0.05L; original: 0 | 光顺系数逐渐恢复至 diffusion 的距离（输入长度单位）；显式设置两参数为 0 时禁用渐增；L 为壁面包围盒最大边长 |
| `--max-skew-angle DEG` | 90 | 瞬时切向/法向推进速度比的限制角度，单位度，0～90；90 不限制，0 纯法向；小于 90 仅支持 tangent 光顺 |
| `--first-layer exact\|march` | height: exact; original: march | 精确法向偏移或方程积分 |
| `--stencil three-point\|compact` | 隐式: three-point；显式: compact | 非均匀周期三点或紧致差分；隐式仅支持三点 |
| `--integrator implicit\|euler\|rk2\|rk4` | implicit | 线性化隐式（默认）或原显式积分器 |
| `--theta T` | 1 | 隐式 theta 参数，0.5 ≤ T ≤ 1；1 为线性化后退 Euler，0.5 为二阶线性化梯形格式 |
| `--cfl C` | 0.25 | 仅显式模式的内部推进步长安全系数，0 < C ≤ 1 |
| `--max-step H` | 0 | 最大内部步长（输入长度单位），0 自动；不限制精确偏移首层 |
| `--max-substeps N` | 1000000 | 单层最大内部子步数；隐式模式限制含失败重试的累计尝试次数 |
| `--format 2d\|3d` | 2d | 2D 或 K=1、Z=0 的 3D Plot3D |
| `--multiblock` | 否 | 添加块数 1 的文件头 |
| `--i-order cw\|ccw` | cw | 输出周向索引，默认 CW 以得到正 Jacobian |
| `--check global\|local` | global | 四角 Jacobian 检查及可选的相邻轮廓全局检查 |
| `--quiet` | 否 | 不输出逐层诊断信息 |
| `--help` | — | 输出帮助 |

### Plot3D 格式

默认二维单块 ASCII：

```text
I J
X(1,1) X(2,1) ... X(I,J)
Y(1,1) Y(2,1) ... Y(I,J)
```

- 输入 M 个独立节点，输出 I=M+1，J=N+1；I 变化最快；X、Y 分量分别连续存储。
- I=1 与 I=M+1 为重复的周期接缝节点；J=1 为壁面；J=N+1 为最外层。
- 内部计算保留逆时针顺序。输出默认反转 I 顺序，但保持起始节点不变，从而满足 `det(∂(x,y)/∂(I,J))>0`。`--i-order ccw` 保留输入顺序，对应负有向 Jacobian。
- `--multiblock` 在维数前写入块数 1。`--format 3d` 写入 I J 1 和全零 Z 分量，适用于需要三维 Plot3D 的读取器。
- 不包含 IBLANK、解变量、Fortran 二进制记录或边界条件文件。

## 3. 数学和离散

记 r=(x,y)，初始壁面弧长参数为 eta，向外推进参数为 xi。原始方程：

```
r_xi = (y_eta, -x_eta) / |r_eta|^(1+alpha) + epsilon * r_etaeta
```

eta 由初始多段线各段长度累计获得，在推进过程中保持不变，包括周期闭合段。程序内部以包围盒最大边长 L 归一化坐标和推进步长，以减轻量纲和尺度差异的影响。`--diffusion` 是无量纲值；相应物理方程中的 epsilon = diffusion * L。

显式积分器默认采用与原 JS 代码相同类型的周期非均匀紧致差分；隐式积分器默认并要求使用非均匀三点差分。紧致差分定义为：

```
betaL f^(d)_(i-1) + f^(d)_i + betaR f^(d)_(i+1)
    = sum(k=-2,...,2) w[k] f_(i+k),  d=1,2
```

用 0 到 6 次 Taylor 项精确匹配确定 7 个未知系数；对均匀点距恢复经典对称六阶紧致格式。对非均匀网格的实际收敛阶数取决于点距分布。周期三对角方程组利用边框消元，预分解一次，每次求解 O(M)，两个坐标分量一起求解。紧致矩阵病态时报告错误，可改用三点格式（非均匀三点二阶导数在任意非光滑点距分布下不保证二阶精度）。

### 推进模式

**height（默认）**

- alpha=0，推进速度的法向分量为 1。
- 默认将光顺投影到切向，即 `epsilon * t * dot(r_etaeta,t)`，避免光顺项直接改变法向推进速度。
- 首层直接使用壁面离散外法向进行精确偏移 `r1=r0+h1*n0`，因此首层节点位移长度等于 h1（浮点舍入误差除外）。法向采用所选差分格式得到的离散切线（隐式默认三点）；尖角没有唯一的连续几何法向。若需首层也按 PDE 积分，设置 `--first-layer march`。
- 后续目标步长 `h_j=h1*q^(j-1)`。因为每层内法向会旋转，切向光顺也改变节点对应关系，最终相对上一层离散法向的节点位移投影仅近似等于 h_j，不保证每个节点的实际层高严格按 q 增长。

**original（原始算法）**

- 默认 alpha=1、完整二阶导数光顺，所有层均通过 PDE 积分。
- 改变 alpha 会改变局部法向速度；特别是 alpha=1 时，无光顺部分满足等参数面积推进：`|det(r_eta,r_xi)|=1`。
- 与 JS 原版的差异：以多段线弦长建立 eta、使用双精度、几何增长替代线性增长、显式模式使用保守自动子步并防止 0 个子步。并非逐点复现原版数据。
- 该模式中 `--first-height` 指第一计算坐标步长，而非严格的物理高度。

例如：

```powershell
.\build\Release\hypergrid.exe --input examples/circle.csv --output original.xyz --first-height 0.002 --growth 1.12 --layers 20 --mode original --diffusion 0
```

### 线性化隐式推进（默认）

设 v=D1 r，z=D2 r，s=|v|，t=v/s，R(x,y)=(y,-x)，P=t t^T，Q=I-P。切向光顺模式的右端项为 F=R t / s^alpha + epsilon P z；完整光顺模式中 P z 替换为 z。

在当前层计算解析 Jacobian A=∂F/∂v，B=∂F/∂z：

```text
A = s^(-1-alpha) R [I-(1+alpha)P]
    + epsilon/s * [(t·z) Q + t (Q z)^T]  (切向光顺)
B = epsilon P
```

完整光顺时 A 不包含第二项，B=epsilon I。用非均匀周期三点差分 D1、D2 构造每个节点的三个 2×2 块（解析线性化包含切向投影的导数，而非冻结投影）：

```text
[I - theta h (A D1 + B D2)] delta_r = h F(r)
r_new = r + delta_r
```

theta=1 为一阶线性化后退 Euler；theta=0.5 为光滑半离散方程的二阶线性化梯形方法。周期块三对角采用边框消元，将最后一个未知块作为两个额外右端项，求解 2×2 Schur 补，再回代。组装及求解均为 O(M)，工作内存 O(M)；全局轮廓检查另计。

默认每个目标层仅求解一次。如果块主元病态、坐标非法、Jacobian 或轮廓检查失败，保持已保存层不变，将该层均分为 2、4、8……个子步重新计算，最多加倍 10 次，且遵守 --max-substeps。所有中间轮廓及最终输出粗层均检查质量。精确偏移首层仍直接计算，不做隐式重试。--max-step 可主动限制子步以控制积分误差；--cfl 只影响显式积分。隐式不是非线性无条件稳定或保正保证。日志 steps 表示成功方案的子步数（精确偏移首层记 1），不包含失败重试。

恢复原显式积分方式可使用 `--integrator rk4`（默认 compact）；如需仅比较时间推进算法，可在显式模式指定 `--stencil three-point`。如需复现旧版不带近壁保护的测试，显式设置 `--smoothing-start 0 --smoothing-end 0`。

### 切向速度限幅（正交性控制）

在 tangent 光顺模式中，令 s=|v|，t=v/s，n=R t，Vn=s^(-alpha)，原始切向速度 Vt0=epsilon*(t·z)。指定角度 phi=--max-skew-angle（单位为度），令 k=tan(phi)：

```text
Vt = clamp(Vt0, -k*Vn, k*Vn)
F = Vn*n + Vt*t
```

仅在超过上限的节点激活限幅；未激活的节点保留原光顺速度。phi=90 表示不限幅（默认，保持旧结果）；phi=0 表示纯法向推进。精确偏移首层不受影响。phi<90 仅支持 tangent 光顺，original 模式需要显式指定 --smoothing tangent 后才能使用。显式及隐式积分器均支持限幅。

对隐式求解，未限幅时使用原解析 Jacobian。限幅激活时令 c=sign(Vt0)*k，F=s^(-alpha)*(R+c*I)*t，因此：

```text
A = s^(-1-alpha) * (R+c*I) * [I-(1+alpha)P]
B = 0
```

phi=0 时取 c=0。切换点除 phi=0 外不可微，代码在等于阈值时选择未限幅分支。组装仍为周期 2×2 块三对角系统，复杂度 O(M)；沿用质量检查和失败时的子步回退。

对前后缘加密的 NACA0012 算例（60 层，h1=1e-4，q=1.15，diffusion=0.002，theta=0.5），建议先对比 5、7.5、10 度：

```powershell
.\build\Release\hypergrid.exe --input examples/naca0012_clustered.csv --output mesh.xyz --first-height 0.0001 --growth 1.15 --layers 60 --theta 0.5 --max-skew-angle 5
```

该参数约束的是半离散方程的瞬时推进方向，而非有限厚度网格单元四角的实际角度。大步长积分、轮廓曲率和周向差分等会引入实际角度偏差。限制过严可能使周向点距分布变差；必要时配合 --max-step 控制时间积分误差。测试数据见 TEST_RESULTS.txt。

### 基于当前点距的切向光顺

设置 `--smoothing-metric current`，将原来的切向光顺速度 epsilon*(t·D2 r) 替换为：

```text
em = r_i - r_(i-1), ep = r_(i+1) - r_i
lm = |em|, lp = |ep|
Vt0 = epsilon * (1/lm - 1/lp)
F = Vn*n + Vt*t
Vt = clamp(Vt0, -Vn*tan(phi), Vn*tan(phi))
```

上述长度和 epsilon 均采用程序内部的归一化量。换算至物理长度时 epsilon_phys = diffusion*L。若两侧长度相等，切向速度为零；否则向较长的一侧移动，趋向均匀点距。点距同比放大时，光顺速度按比例减小。不直接采用当前弧长下的二阶导数，因为精确弧长参数的二阶导数为法向曲率，无法有效驱动切向重分布。

法向项及 D1（固定初始壁面弧长坐标）不变，精确偏移首层不变。每个隐式子步和显式 RK 阶段使用当前坐标重新计算实际点距，不做跨层冻结。默认 wall 保留原光顺公式；若需恢复未启用近壁保护的旧结果，同时设置 --smoothing-start 0 --smoothing-end 0。current 仅适用于 tangent 模式，不支持 full。diffusion=0 或 max-skew-angle=0 时仍为纯法向推进。

未限幅时，令 gm=∂Vt0/∂r_(i-1)，gp=∂Vt0/∂r_(i+1)，则：

```text
gm = epsilon*em/lm^3
gp = epsilon*ep/lp^3
g0 = -gm-gp
A = s^(-1-alpha)*R*[I-(1+alpha)P] + (Vt0/s)*Q
Gm = t*gm^T
Gp = t*gp^T
```

设 D1 的三点系数为 a-, a0, a+，新右端项对相邻坐标的 Jacobian 块为 A*a-+Gm、A*a0-Gm-Gp、A*a++Gp。若限幅激活，沿用此前限幅分支的 A，Gm=Gp=0。组装及求解仍为 O(M) 的周期 2×2 块三对角系统；显式积分器相应采用当前点距梯度和切向方向导数的局部范数估计调整子步（启发式稳定性控制）。

推荐先尝试加密 NACA0012 60 层算例：

```powershell
.\build\Release\hypergrid.exe --input examples/naca0012_clustered.csv --output mesh.xyz --first-height 0.0001 --growth 1.15 --layers 60 --theta 0.5 --smoothing-metric current --diffusion 0.03
```

current 与 wall 的光顺强度不能直接用相同 diffusion 值比较。该算例中 current 的 0.02～0.03 较为均衡，外层平均和最大正交偏差及外层周向点距比优于旧算法限幅 5 度。该结论仅针对测试算例。过大的 diffusion 会削弱原有近壁加密并使锐尾缘单元质量下降；0.2 且不限幅的测试在第 2 层失败。严格限幅和强光顺同时使用可能导致更多回退，并不保证更好的点距分布。不要将瞬时速度限幅解释为单元角度的硬约束。完整测试见 TEST_RESULTS.txt。

### 近壁保护：按名义推进距离逐步恢复光顺

新 current 光顺的倒数点距差在加密壁面附近可能很大，不宜直接在全域使用较强的 diffusion。height 模式默认开启近壁保护：如果未指定任一距离参数，自动设置 a=0.005L、b=0.05L，L 为壁面包围盒最大边长。original 模式默认 a=b=0，以保持原算法兼容性。显式指定任一距离参数时，按原命令行语义使用给定值，另一参数未指定则为 0；要求 0<=a<b 或两者均为 0（关闭保护）。设 H 为累计名义推进距离，令

```text
u = clamp((H-a)/(b-a), 0, 1)
epsilon_eff = diffusion * u^2 * (3-2*u)
```

在每个输出层的名义推进距离中点取 H，并在该层的隐式子步及显式 RK 各阶段保持该预定系数不变。实际点距、法向和其解析 Jacobian 仍按原算法随坐标更新。对 wall/current 和所有支持的积分器均有效，保留角度限幅及质量检查。逐层日志新增 effective-diffusion。距离参数使用输入坐标长度单位，height 模式下是名义距离而非逐点真实壁面距离；original 模式中仍为推进参数。渐增的系数冻结在输出层尺度，max-step 细分不会进一步细分这个系数。首层精确偏移不受影响。

对于单位弦长加密 NACA0012，60 层、h1=1e-4、q=1.15，默认自动距离即为 0.005/0.05，无需再显式指定，推荐：

```powershell
.\build\Release\hypergrid.exe --input examples/naca0012_clustered.csv --output mesh.xyz --first-height 0.0001 --growth 1.15 --layers 60 --theta 0.5 --smoothing-metric current --diffusion 0.03
```

前 15 层坐标与纯法向推进完全相同，前 10 层平均四角正交偏差从 10.7265 度降至 1.1250 度，外侧 10 层平均/最大偏差 2.298/8.706 度，最外层周向点距比 18.046，60 层仍只需 60 次推进（含首层偏移）。自动距离是几何尺度启发式默认值，并非物理边界层厚度模型；参数需结合实际边界层厚度选择。首层仍沿离散节点法向精确偏移；新增日志 Wall orthogonality (deg)，分别报告与离散节点切线的最大偏差 nodal max，以及与相邻壁面线段的平均和最大偏差 face mean/max。smoothing-ramp=auto/manual/off 表示实际采用的控制方式。该算例 nodal max 约为 3.12e-11 度，face mean/max 为 1.125/81.7385 度。注意单块 O 型锐尾缘共用节点的两侧壁面切向并不共线，通常不可能使同一条法向网格线同时与两侧壁面正交；该改动不能消除锐尾缘固有的局部角度误差。

### 稳定性与质量

显式推进步长结合一阶/二阶离散算子的保守无穷范数上界与局部切线模长估计；推荐 RK4。显式 Euler/RK2 对中心差分的某些波动模态不具有无条件稳定性，不能仅靠减小 CFL 保证所有复杂网格都稳定。此处步长控制是保守启发式而不是非线性严格稳定性证明。

每层输出前检查所有四边形的四个角点的归一化 Jacobian（-cross/(两边长度乘积)，以内部 CCW 索引的外向方向为正），并检查节点法向位移为正。默认同时检查外轮廓自交、相邻轮廓相交及嵌套。逐层报告目标步长、实际法向投影的最小/平均/最大值、最小归一化 Jacobian 和子步数。检查失败且隐式自动重试仍不能生成合格层（或显式检查失败）时，返回非零退出码，指明层数及问题位置，在网格全部生成前不会写入输出文件；不执行网格修补或自动截断。若目标输出文件已经存在，失败时原文件不会被删除。

强烈聚集的点距使显式扩散步长很小；可调整 `--diffusion`、重新分布壁面点或选择更稳健的三点离散。不能保证任意深凹轮廓、极尖尾缘、任意远场距离均能成功生成正 Jacobian 网格。

## 4. 示例与回归测试

`examples/circle.csv`：120 个独立圆周点（半径 0.5）。
`examples/naca0012.csv`：160 个独立点，闭合锐尾缘 NACA0012，按多段线弧长近似等距重采样。
`examples/naca0012_clustered.csv`：160 个独立节点，前后缘加密（90% 余弦、10% 等弧长混合）。
`examples/naca0012.xyz`：原显式 RK4/compact 版本的均匀翼型 40 层参考网格（保留用于对比）。

```powershell
ctest --test-dir build -C Release --output-on-failure
```

本次使用 Linux GCC/CMake Release 构建并通过自动测试；未在 Windows/MSVC 或 Tecplot 中实测。测试包含：圆形等速法向推进解析解、alpha=1 无光顺圆形解析解、NACA0012 首层高度与全网格四角 Jacobian、非均匀周期节点、替代积分器和三点差分、三维单层多块格式、周期接缝、输入节点保留及非法输入拒绝。

测试结果见 `TEST_RESULTS.txt`，其中追加了隐式算法的收敛阶验证及加密 NACA0012 40/60 层性能测试。60 层默认隐式只需 59 次线性求解和一次首层偏移；相对原 RK4/compact 实测约 90 倍加速。大步长的一阶时间积分会改变外层网格位置，如需更接近显式小步长解，建议 `--theta 0.5` 或使用 `--max-step` 控制步长。NACA0012 算例虽未出现负 Jacobian，锐尾缘附近最小归一化 Jacobian 约 0.143，说明仍存在较大局部倾斜；并非已针对所有翼型完成网格质量优化。

## 5. 文件组织

```text
CMakeLists.txt
src/main.cpp           命令行、读写、差分、显式推进、质量检查
src/implicit.hpp       解析 Jacobian、周期 2×2 块三对角、隐式推进及回退
examples/*.csv        测试输入
examples/naca0012.xyz 生成的示例网格
tests/regression.py   自动回归测试
TEST_RESULTS.txt      实测结果
README.md             使用和算法说明
```

移植后续若需集成到 CFD 求解器，可将 `Derivative`、`Marcher` 与 `Ring/Grid` 抽离为独立库，并将输入和 Plot3D 输出保留在命令行前端。
