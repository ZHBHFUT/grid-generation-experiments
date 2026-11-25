# HyperGrid：二维双曲型 O 型结构网格生成器

独立 C++17 命令行程序。基于 Ricky Reusser 的 [grid-generation-experiments](https://github.com/rreusser/grid-generation-experiments) 中的双曲推进方程、周期紧致差分及显式积分算法重新实现，输入由 NACA 参数改为任意简单闭合的逆时针二维多段线。原始算法说明：[Hyperbolic Grid Generation](https://rreusser.github.io/2016/02/07/hyperbolic-grid-generation/)。上游 `package.json` 声明 ISC 许可证；上游仓库未包含单独的 LICENSE 文件。本项目未复制原网页、npm 依赖或第三方库。

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
| `--first-layer exact\|march` | height: exact; original: march | 精确法向偏移或方程积分 |
| `--stencil compact\|three-point` | compact | 周期非均匀紧致差分或显式三点差分 |
| `--integrator euler\|rk2\|rk4` | rk4 | 显式积分方式，推荐 RK4 |
| `--cfl C` | 0.25 | 内部推进步长安全系数，0 < C ≤ 1 |
| `--max-step H` | 0 | 最大内部步长（输入长度单位），0 自动 |
| `--max-substeps N` | 1000000 | 单层最大内部子步数 |
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

一阶和二阶导数采用与原 JS 代码相同类型的周期非均匀紧致差分：

```
betaL f^(d)_(i-1) + f^(d)_i + betaR f^(d)_(i+1)
    = sum(k=-2,...,2) w[k] f_(i+k),  d=1,2
```

用 0 到 6 次 Taylor 项精确匹配确定 7 个未知系数；对均匀点距恢复经典对称六阶紧致格式。对非均匀网格的实际收敛阶数取决于点距分布。周期三对角方程组利用边框消元，预分解一次，每次求解 O(M)，两个坐标分量一起求解。紧致矩阵病态时报告错误，可改用三点格式（非均匀三点二阶导数在任意非光滑点距分布下不保证二阶精度）。

### 推进模式

**height（默认）**

- alpha=0，推进速度的法向分量为 1。
- 默认将光顺投影到切向，即 `epsilon * t * dot(r_etaeta,t)`，避免光顺项直接改变法向推进速度。
- 首层直接使用壁面离散外法向进行精确偏移 `r1=r0+h1*n0`，因此首层节点位移长度等于 h1（浮点舍入误差除外）。法向采用紧致差分或所选三点差分得到的离散切线；尖角没有唯一的连续几何法向。若需首层也按 PDE 积分，设置 `--first-layer march`。
- 后续目标步长 `h_j=h1*q^(j-1)`。因为每层内法向会旋转，切向光顺也改变节点对应关系，最终相对上一层离散法向的节点位移投影仅近似等于 h_j，不保证每个节点的实际层高严格按 q 增长。

**original（原始算法）**

- 默认 alpha=1、完整二阶导数光顺，所有层均通过 PDE 积分。
- 改变 alpha 会改变局部法向速度；特别是 alpha=1 时，无光顺部分满足等参数面积推进：`|det(r_eta,r_xi)|=1`。
- 与 JS 原版的差异：以多段线弦长建立 eta、使用双精度、几何增长替代线性增长、使用保守自动子步并防止 0 个子步。并非逐点复现原版数据。
- 该模式中 `--first-height` 指第一计算坐标步长，而非严格的物理高度。

例如：

```powershell
.\build\Release\hypergrid.exe --input examples/circle.csv --output original.xyz --first-height 0.002 --growth 1.12 --layers 20 --mode original --diffusion 0
```

### 稳定性与质量

推进步长结合一阶/二阶离散算子的保守无穷范数上界与局部切线模长估计；默认 RK4。显式 Euler/RK2 对中心差分的某些波动模态不具有无条件稳定性，不能仅靠减小 CFL 保证所有复杂网格都稳定。此处步长控制是保守启发式而不是非线性严格稳定性证明。

每层输出前检查所有四边形的四个角点的归一化 Jacobian（-cross/(两边长度乘积)，以内部 CCW 索引的外向方向为正），并检查节点法向位移为正。默认同时检查外轮廓自交、相邻轮廓相交及嵌套。逐层报告目标步长、实际法向投影的最小/平均/最大值、最小归一化 Jacobian 和子步数。检查失败立即返回非零退出码，指明层数及问题位置，在网格全部生成前不会写入输出文件；不执行网格修补或自动截断。若目标输出文件已经存在，失败时原文件不会被删除。

强烈聚集的点距使显式扩散步长很小；可调整 `--diffusion`、重新分布壁面点或选择更稳健的三点离散。不能保证任意深凹轮廓、极尖尾缘、任意远场距离均能成功生成正 Jacobian 网格。

## 4. 示例与回归测试

`examples/circle.csv`：120 个独立圆周点（半径 0.5）。
`examples/naca0012.csv`：160 个独立点，闭合锐尾缘 NACA0012，按多段线弧长近似等距重采样。
`examples/naca0012.xyz`：上述翼型，h1=1e-4，q=1.15，40 层，默认 height 模式。

```powershell
ctest --test-dir build -C Release --output-on-failure
```

本次使用 Linux GCC/CMake Release 构建并通过自动测试；未在 Windows/MSVC 或 Tecplot 中实测。测试包含：圆形等速法向推进解析解、alpha=1 无光顺圆形解析解、NACA0012 首层高度与全网格四角 Jacobian、非均匀周期节点、替代积分器和三点差分、三维单层多块格式、周期接缝、输入节点保留及非法输入拒绝。

测试结果见 `TEST_RESULTS.txt`。NACA0012 算例虽未出现负 Jacobian，锐尾缘附近最小归一化 Jacobian 约 0.143，说明仍存在较大局部倾斜；并非已针对所有翼型完成网格质量优化。

## 5. 文件组织

```text
CMakeLists.txt
src/main.cpp           命令行、读写、紧致差分、周期三对角、推进、质量检查
examples/*.csv        测试输入
examples/naca0012.xyz 生成的示例网格
tests/regression.py   自动回归测试
TEST_RESULTS.txt      实测结果
README.md             使用和算法说明
```

移植后续若需集成到 CFD 求解器，可将 `Derivative`、`Marcher` 与 `Ring/Grid` 抽离为独立库，并将输入和 Plot3D 输出保留在命令行前端。
