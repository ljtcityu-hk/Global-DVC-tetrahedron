# 从用户网格文件求解 Tet4 Global DVC

本目录独立运行，不依赖上一级 `tetrahedron` 目录。主程序和 demo 只读取用户已有网格，绝不生成网格或重新进行 Delaunay 剖分。
支持任意非规则节点位置和有效共形四节点线性四面体（Tet4）连接关系，包括凹域和内部孔洞。
“任意”指节点坐标及网格拓扑；当前文件解析支持下面明确列出的 MAT、CSV/TXT 格式，不代表自动识别所有网格软件格式。Tet10、混合单元、Abaqus INP/Gmsh 等原生文件需先转换。

## 直接运行样例

在 MATLAB 切换到本目录：

```matlab
demo_userMesh_tetrahedron
```

Section 1 设置两份图像 MAT 和一份已有网格 MAT 路径。切换注释中的 `meshSource` 可改用节点/单元 CSV 文件。
Section 2 设置 alpha、初值和坐标转换；无需 ROI 节点网格或 winstepsize。
运行后保存 `results_userMesh_tetrahedron.mat` 和 `displacement_userMesh_tetrahedron.csv`。
CSV 四列为原始节点 ID、u、v、w。结果中的节点顺序与输入节点表一致，单元结果顺序与输入单元表一致。
依赖 MATLAB 和 Image Processing Toolbox，绘图无需 PDE Toolbox。

## 网格文件规范

MAT 文件变量（或同名字段放在一个 `DVCmesh` 结构体中）：

```matlab
coordinatesFEM  % N×3 double，各行 [x,y,z]
elementsFEM     % M×4，每行四个节点编号
nodeIDs         % 可选 N×1，非连续编号也可；省略则为 1:N
elementIDs      % 可选 M×1，省略则为 1:M
save('mesh.mat','coordinatesFEM','elementsFEM','nodeIDs','elementIDs');
```

**存在 nodeIDs 时，elementsFEM 引用这些 ID；没有 nodeIDs 时引用 MATLAB 的 1-based 节点行号。**
仅有节点坐标无法确定已有四面体及孔洞，必须同时提供单元连接关系。

CSV/TXT 使用两份纯数值文件，无表头、无计数行：

- 节点文件：三列 `x,y,z`，或四列 `nodeID,x,y,z`。
- 单元文件：四列 `n1,n2,n3,n4`，或五列 `elementID,n1,n2,n3,n4`。
- 逗号或空白分隔。无 ID 时节点编号为 1:N；显式 ID 可以从 0 开始，也可以不连续。

```matlab
meshSource = struct('nodesFile','nodes.csv','elementsFile','elements.csv');
result = main_FE_GlobalDVC_userMesh_tetrahedron('ref.mat','def.mat',meshSource,DVCpara,U0);
```

读取器检查未知 ID、重复 ID/坐标/节点/单元、退化单元、未使用节点及非流形共享面，并修正负方向四面体（交换局部第 2、3 节点）。输入单元 ID 和节点行顺序保留。
输入网格须有效、共形、无几何自交；当前检查不包含任意两个不共享面的四面体之间的完整几何相交检测。悬挂节点或裂缝重合节点不在本阶段支持范围。

## 图像与坐标

图像 MAT 与原 GenerateVolMatfile 一致，含 `vol{1}` 三维灰度数组；也接受直接存成数值三维数组的 `vol`。
`x,y,z` 对应 MATLAB 图像第一、第二、第三维。已转换的 vol 不应再次 permute。

默认节点坐标为 1-based 图像体素坐标。物理坐标转换为：

```matlab
DVCpara.meshScale  = 1 ./ voxelSpacing;
DVCpara.meshOffset = 1 - imageOrigin ./ voxelSpacing;
% X_voxel = X_input .* meshScale + meshOffset
```

此处 imageOrigin 是首个体素中心的物理坐标，轴方向必须已对齐。旋转坐标系需用户先转换节点。
网格所有节点应在 `[4,size(Img,axis)-3]` 的七点梯度有效域中，位移后采样位置也须在变形图像内。
U0 与输出 U 均为体素位移，按输入节点行顺序 `[u1,v1,w1,u2,v2,w2,...]'` 排列；省略或传 `[]` 则用零初值。
物理位移可由每个方向的体素位移乘 voxelSpacing 得到。非等距体素下输出的梯度/应变仍基于体素坐标，物理应变应先按尺度变换梯度后再计算。

## 求解与结果

使用输入连接关系给整数体素分配唯一单元，只采样网格内部。归一化统计量同样排除孔洞及域外体素。
灰度项按单位体素求和，正则项为 `alpha * integral |grad(u)|^2 dV`，按四面体体积积分。
细到没有分配体素的单元会报告警告和 `info.elementSampleCounts`，其结果依赖邻域及正则约束；弱纹理或欠约束的矩阵会明确报错。

`F` 为节点位移梯度（不是 I+grad(u)），九分量按列排列 `[ux,vx,wx,uy,vy,wy,uz,vz,wz]`；`Strain` 为对应对称应变。
`StrainType=0/1/2`：小应变/Euler-Almansi/Green-Lagrange。
`elementGradient`、`elementStrain` 每行对应一个输入单元。节点梯度通过体积加权恢复。
`DVCmesh.inputNodeIDs/inputElementIDs` 保留输入编号，内部 `elementsFEM` 已转换为节点行号，不可误作原始 ID。
`info.converged/status/objective` 用于确认收敛；本阶段一次计算一对图像，无外加位移边界条件或多帧增量流程。

## 样例及验证

`sample_userMesh_tetrahedron` 包含同一网格的 MAT 和 CSV 版本：非整数扰动节点、非连续编号、打乱节点行顺序和部分反向单元，并保留中心孔洞。
图像为 40×40×40 uint16 合成体图像，真值平移 `[0.35,-0.25,0.20]` 体素。
网格文件已经生成并随目录提供，demo 不需要也不调用任何网格生成器。
运行 `report = test_userMesh_tetrahedron` 检查两种文件格式、编号映射、方向修正、孔洞、仿射梯度、异常输入及端到端位移恢复。

本机 MATLAB 已通过上述测试：64 个节点、156 个 Tet4，13041 个唯一体素，78 个反向单元得到修正。
不归一化的合成平移回归测试 RMSE 约 0.0026 体素；默认 demo 启用独立图像归一化，RMSE 约 0.0200 体素，平均位移约 `[0.3476,-0.2558,0.2008]`。
两者均从零初值收敛。孔洞排除、绘图、MAT 保存及带原节点 ID 的 CSV 导出均已验证。
归一化在有限区域内会改变两幅图像的统计量，因此其误差与不归一化测试不同；实际数据应按亮度变化和纹理选择参数。
以上为合成样例验证，用户自己的真实 CT 和网格仍需单独验证。
