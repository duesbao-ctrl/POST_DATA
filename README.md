# POST_DATA（非界面版）

面向 LAMMPS chunk、cluster、mass-v、mass-x 和二维孔隙网络数据的 MATLAB
后处理核心。该项目只提供脚本/API，不包含图形界面；计算、输出和绘图能力与
`POST_DATA2` 共用同一套模块化实现，兼容 MATLAB R2016b，不使用 `+package`。

Windows 7 / MATLAB R2016b 下的结果视图中文文本由
`src/core/resources/ui_zh_CN.tsv` 以 UTF-8 显式读取；可执行 `.m` 文件
保持纯 ASCII，不依赖 Windows 系统代码页。

## 启动与基本用法

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA')
postdata_startup

request = pd_create_request('mass-x');
request.filePath = fullfile(pwd, 'fixtures', 'generated', ...
    'large_bin1d_dx_0.025.txt');
request.analysisOptions = {'ChunkDim', '1d', 'SphDimension', 2, ...
    'InitialDensity', 7.3, 'ParticleSpacing', 0.005, ...
    'RawLengthUnitUm', 10, 'TransverseWidth', 2.46, ...
    'CoordinateFactor', 10};
request.makePlots = true;
request.plotOptions = struct('FontName', 'Times New Roman', ...
    'LineWidth', 1.8, 'SeriesStyleMode', 'color-and-style');
result = postdata_run(request);
```

`postdata_run` 是唯一计算调度入口。旧脚本可继续使用 `run_analysis`，该函数现在
只是新请求 API 的兼容包装器，不再保存第二套分析实现：

```matlab
result = run_analysis('chunk', ...
    'BaseDir', fullfile(pwd, 'fixtures', 'sample'), ...
    'ChunkFile', 'bin1d_strain_rate_dx_1.txt', ...
    'ChunkDim', 'auto', 'Variable', 'strainRate', ...
    'DoPlot', false, 'ProgressMode', 'off');
```

## 分析能力

- `chunk`：自动识别一维/二维；支持原始场、温度/速度/应力等派生场、任意一维
  场梯度以及 `du/dx + (u/rho)*d(rho)/dx` 应变率。
- `cluster`：团簇直径、直方图/PDF/CDF、拟合、位置分箱统计；支持简单立方、
  显式三维粒子体积和薄层准二维投影三种物理尺寸模型。
- `mass-v`（规范类型 `vx`）：速度—面密度微分/累积分布，可修改速度列、换算
  系数、单位、密度列和累积方向。
- `mass-x`（规范类型 `massx`）：仅处理 SPH chunk 的 `Ncount`。根据初始密度、
  粒子间距和横向宽度推导局部/累计面密度，从大 x 向小 x 累加，单位
  `mg/cm^2`；支持 bin1d、bin2d 全宽和多个 y 分片。
- `network2d`：孔隙/基体几何、连通分量、数字拓扑、方向剖面、演化、PLIC
  切割单元和可选骨架形态统计。

`mass-v`、`mass_v`、`mass-vx`、`massvx` 均映射到 `vx`；`mass-x` 与
`mass_x` 映射到 `massx`。

## 新版 SPID chunk 兼容

读取层同时兼容原有三行表头与新版 SPID 可扩展表头。新版文件中的
`# Units`、`# Name ... Kind ...` 和逐帧 `# Time` 会写入预检及结果元数据。

- `SelectBy=Time` 优先使用文件内的物理时间，不再强制依赖 Slurm；旧文件仍按
  原逻辑使用 Slurm 映射。
- `spatial` 自动识别一维、二维和三维坐标。
- `field` 类型的 mass-v 自动使用 `Coord1`；旧格式仍优先使用 `Chunk`。
- `cluster` 同时支持新版 `x/y/z` 和旧版 `c_x/c_y/c_z`。
- `UseFileUnitMetadata=true` 时，已知 SPID 单位制会自动换算 mass-v 和 mass-x；
  没有单位元数据的旧文件继续使用原有手动换算参数。

## SPH mass-x 单位与分片

`ChunkDim` 是 chunk 的一维/二维分箱类型，`SphDimension` 是 SPH 模拟本身的
二维/三维物理维度。二维 SPH 使用 `rho0*d0^2` 的单位厚度粒子质量并除以 y
统计宽度；三维 SPH 使用 `rho0*d0^3` 的粒子质量并除以 `Wy*Wz`，因此还需
设置 `OutOfPlaneWidth`。

`InitialDensity` 的单位为 `g/cm^3`。`ParticleSpacing`、`TransverseWidth`、
`OutOfPlaneWidth`、`SliceCentersY` 和 `SliceWidthsY` 均使用原始 SPH 坐标
单位。通常一个原始坐标单位等于 `10 um`，应设置 `RawLengthUnitUm=10`；
绘图坐标另设 `CoordinateFactor=10` 后显示为 `um`。这两个参数分别控制物理
换算和显示，不能混用。

bin1d 必须指定完整 y 宽度。bin2d 默认使用完整 y 网格，也可通过
`SliceCentersY` 与 `SliceWidthsY` 同时输出多个位置/宽度的曲线。边界落在
y 网格内部时，按网格覆盖比例分配 `Ncount`。局部面密度允许随 x 波动；从大
x 向小 x 求和后的累计面密度按 x 升序显示时严格单调不增。

SPH 原始场和 SPH mass-x 不使用 `dV`。`dV` 参数只属于 MD 守恒量的派生分析
路径，各 MD 派生公式按自身需要使用该值。

## 论文级绘图与导出

所有 `postdata_run` 独立图和 `pd_render_result` 调用共享
`pd_plot_option_catalog`。可修改字体、字号、标题、坐标标签、LaTeX/TeX 解释器、
线宽/线型/标记、色盲配色、自定义 RGB 色序、彩色+线型或纯黑白线型、网格、
图例、色条、色限、坐标范围、对数轴、反向轴和等比例轴。

```matlab
request.plotOptions = struct( ...
    'SeriesStyleMode', 'monochrome', ...
    'TextInterpreter', 'latex', ...
    'FontSize', 11, 'TitleFontSize', 12, ...
    'LineWidth', 1.8, 'ShowGrid', false);
```

`pd_output_option_catalog` 支持 MAT、摘要/明细 CSV、PNG、FIG、PDF 和可复现清单；
默认 PNG 为 300 DPI，`FigureWidthCm` 与 `FigureHeightCm` 可直接按期刊栏宽设置。
默认输出集中到项目的 `outputs/`，多格式导出会先完整写入隔离临时目录，任一格式
失败时不留下半成品。

每幅结果图的实际数值由 `pd_result_plot_data` 统一生成，直方图、CDF、方向均值、
mass-v/mass-x 曲线和二维网格均可直接查看或导出：

```matlab
viewData = pd_result_plot_data(result, 'cumulative');
disp(viewData.ColumnNames)
disp(viewData.Values)
pd_write_plot_data_csv('mass_x_cumulative.csv', viewData);
```

输出选项 `SavePlotDataCSV=true` 会为当前结果的每个可用视图分别生成一份 CSV，
其数据与 `pd_render_result` 使用的数值完全一致。

## 稳定性与诊断

- 脚本请求与界面版共用参数目录验证，未知、重复或非法计算/绘图参数在文件读取前报错。
- 输入预检验证首块结构、数据列数和归一化后的重复字段名。
- 长文件索引支持取消；损坏/过期缓存自动重建，缓存采用原子替换。
- 网络演化使用内存索引，不要求输入目录可写。
- 运行日志超过 5 MB 自动轮换，可通过 `pd_system_diagnostics()` 检查安装与可选工具箱。

```matlab
report = pd_system_diagnostics();
disp(report)
```

## 目录结构

```text
POST_DATA/
  postdata_run.m              统一请求调度入口
  run_analysis.m              旧调用方式兼容包装器
  postdata_startup.m          路径初始化
  src/
    core/                     请求、参数目录、验证、日志、UTF-8 文本资源
    io/                       chunk/Slurm 读取、索引、输入预检
    analysis/                 数值分析与 network2d 子模块
    plot/                     逐视图数值契约、论文样式和多视图渲染
    export/                   原始明细、逐视图数据和图形导出
  examples/                   快速开始和大数据验证
  fixtures/sample/            小型确定性回归数据
  fixtures/generated/         一维、二维、颗粒、mass-v 大型数据
  tests/test_postdata.m        自动化测试
```

项目中不存在 `postdata_app.m` 或 `src/app`；需要界面时启动同级 `POST_DATA2`。

## 测试

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA')
postdata_startup
results = runtests('tests');
assert(numel(results) > 0, 'No tests were discovered.');
assertSuccess(results);
```

GUI 专属的 5 个测试在本项目中会被明确标记为跳过，其余计算、读取、绘图、导出、
大数据 `mass-v/mass-x`、梯度/应变率和物理颗粒尺寸测试与 `POST_DATA2` 相同。

大数据验证：

```matlab
addpath('examples')
paths = generate_large_test_data();
report = run_large_verification();
```

默认验证图和摘要写入 `outputs/verification/`，该目录是可随时
重建的运行输出，不纳入源码和分发包。

## 同步约束

`POST_DATA2` 是包含 GUI 的上游超集。非界面白名单通过
`POST_DATA2/tools/sync_non_gui_to_post_data.ps1` 同步；该脚本显式排除界面目录，
删除根目录旧重复实现和 MATLAB 临时缓存，避免两套计算代码再次漂移。
