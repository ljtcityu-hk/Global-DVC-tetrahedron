# Tetrahedral Global Digital Volume Correlation

基于用户输入四面体网格的三维全局数字体相关（Global DVC）MATLAB 程序。

## 项目用途声明

本项目用于从参考状态与变形状态的三维灰度图像中估计材料内部的三维位移场及应变场。程序采用四节点线性四面体（Tet4）作为位移插值单元，通过全局灰度匹配和位移梯度正则化求解节点位移，再计算单元位移梯度、恢复节点梯度并输出应变结果。

本项目在 Jin Yang 的 FE_Global_DVC 代码基础上进行扩展，主要面向用户已经建立的四面体体网格。求解入口直接读取节点坐标和单元连接关系，使计算区域能够描述非规则形状、凹域及内部孔洞，并保持相邻单元之间的位移连续性。

代码可用于数字体相关算法研究、材料变形测量方法开发、合成体图像验证，以及与外部网格处理流程的集成。对于 CT 等实际体图像，应根据图像纹理、噪声、体素尺寸和预期位移选择网格与求解参数，并通过相应实验或已知变形数据验证结果。

这里的有限元网格用于参数化位移场。当前程序不包含材料本构模型、外力加载或应力求解，输出属于基于图像相关的运动学测量结果。

**English overview:** This MATLAB project estimates three-dimensional displacement and strain fields from a pair of volumetric grayscale images using finite-element-based global digital volume correlation. It extends Jin Yang's FE_Global_DVC workflow to user-supplied, conforming four-node tetrahedral meshes. The solver reads existing node coordinates and element connectivity, minimizes an image-matching objective with displacement-gradient regularization, and exports nodal displacements and recovered deformation measures. It is intended to support DVC research, validation, and integration with external mesh workflows.

## 主要功能

- 读取 MAT 体图像，以及 MAT 或 CSV/TXT 格式的 Tet4 网格。
- 支持非规则节点位置、非连续节点编号，并保留输出与输入编号的对应关系。
- 检查常见网格输入错误，并修正四面体的负方向节点排序。
- 按输入网格确定体素归属，避免共享面重复采样，保留孔洞及实际域边界。
- 求解全局节点位移，输出单元梯度、体积加权节点梯度，以及小应变、Euler–Almansi 或 Green–Lagrange 应变。
- 提供网格表面和剖开结果显示、MAT 结果保存及节点位移 CSV 导出。

## 快速开始

需要 MATLAB 和 Image Processing Toolbox。进入本目录后，选择一个示例运行：

```matlab
demo_userMesh_tetrahedron           % 带孔洞的非规则节点网格
demo_large_userMesh_tetrahedron     % 128³ 图像：平移和仿射变形
demo_cylinder_userMesh_tetrahedron  % 圆柱体：平移和轴向拉伸
```

样例所需的图像和网格文件已提供。示例读取已有网格；各样例目录中的生成程序仅用于离线复现数据，不在 DVC 求解过程中调用。

使用自己的数据时，修改 demo 的图像及网格文件路径，并按实际坐标单位设置参数。输入格式、编号规则、坐标转换及输出字段见[详细使用说明](README_userMesh_tetrahedron.md)。

## 当前范围与验证情况

当前版本支持有效共形的四节点线性四面体体网格，一次处理一对图像。节点文件必须配有单元连接关系；Tet10、混合单元及 INP/Gmsh 等其他原生文件格式尚未直接接入。

默认初始位移为零，也可由用户传入节点初值。目前未集成 FFT 初始位移搜索、多帧增量流程或外加位移边界条件。较大位移需要合适的初值，不能仅根据小位移样例的收敛情况推断其可解范围。

图像坐标采用 `Img(x,y,z)`，位移默认以体素为单位。非等距体素或外部物理坐标应按详细说明处理，尤其要区分体素坐标下的梯度与物理坐标下的梯度。

已完成网格连通性、孔洞保留、编号映射、仿射梯度和合成图像位移恢复等检查。提供的误差与耗时来自合成纹理样例，不能直接视为真实 CT 数据的精度或性能保证。

- [基础回归测试](test_userMesh_tetrahedron.m)
- [128³ 大算例说明](sample_large_userMesh_tetrahedron/README_large_tetrahedron.md)
- [圆柱算例说明](sample_cylinder_tetrahedron/README_cylinder_tetrahedron.md)

## 来源与许可证

原始 FE_Global_DVC 代码版权归 Jin Yang；本项目修改及新增内容的版权声明为 Jiatai。来源和修改范围见 [NOTICE.md](NOTICE.md)。

```text
Copyright (c) 2020, Jin Yang.
Copyright (c) 2026, Jiatai.
```

本项目按 [BSD 2-Clause License](license.txt) 发布。复制、修改及再分发时，应遵守该许可证并保留相应版权声明、许可条件及免责声明。本 README 的用途说明不增加许可证之外的用途限制。
