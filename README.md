# POST_DATA2

面向 LAMMPS chunk、cluster、mass-v、mass-x 和二维孔隙网络数据的 MATLAB 后处理软件。项目仅保留一套模块化实现，兼容 MATLAB R2016b，不依赖 App Designer，也不使用 `+package` 目录。

Windows 7 / MATLAB R2016b 的中文界面由 `src/core/resources/ui_zh_CN.tsv` 以 UTF-8 显式加载；可执行 `.m` 文件保持纯 ASCII，不依赖系统代码页，因此不会因英文或中文区域设置产生乱码。

## 启动

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA2')
postdata_app
```

脚本方式需要先初始化路径：

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA2')
postdata_startup

request = pd_create_request('massx');
request.filePath = fullfile(pwd, 'fixtures', 'sample', 'bin1d_dx_0.5.txt');
request.analysisOptions = {'ChunkDim', '1d', ...
    'SphDimension', 2, ...
    'InitialDensity', 7.3, ...       % g/cm^3
    'ParticleSpacing', 0.1, ...      % raw SPH coordinate units
    'RawLengthUnitUm', 10, ...       % 1 raw unit = 10 um
    'TransverseWidth', 1, ...        % full y width, raw units
    'CoordinateFactor', 10};         % plot x in um
request.makePlots = false;
result = postdata_run(request);
```

## 分析类型

- `chunk`：一维或二维原始场/派生场。`ChunkDim=auto` 会根据 `Coord2/c_y` 是否存在自动识别维度，避免一维文件误报缺少 `Coord2`。
- `cluster`：团簇尺寸、概率分布、CDF、拟合和位置分箱统计。
- `vx`：mass-v 面密度分布，默认从最大速度向最小速度累积。
- `massx`：仅用于平面二维 SPH 结果，根据 chunk 中的 `Ncount`、初始密度、初始粒子间距和 y 统计宽度计算 mass-x；固定从最大 x 向最小 x 累积，纵轴单位为 `mg/cm^2`。
- `network2d`：二维孔隙/基体几何、拓扑、连通性、方向剖面和演化统计。

### SPH mass-x 物理模型

SPH 的 `bin1d` / `bin2d` 文件不需要、也不应包含面密度列。程序读取每个
chunk 的 `Ncount`。`ChunkDim` 表示 chunk 分箱维度，`SphDimension` 表示
模拟本身是二维还是三维，两者不能混为一谈。对于平面二维 SPH，单粒子代表的
单位厚度质量为

```text
particleLineMass = rho0 * (d0 * Lunit)^2              [g/cm]
localArealDensity = Ncount * particleLineMass
                    / (Wy * Lunit) * 1000             [mg/cm^2]
```

其中 `rho0=InitialDensity`，单位为 `g/cm^3`；`d0=ParticleSpacing` 和
`Wy` 都使用原始模拟坐标单位；`Lunit=RawLengthUnitUm*1e-4 cm`。本 SPH
单位制通常为一个原始坐标单位等于 `10 um`，所以 `RawLengthUnitUm=10`，
同时绘图坐标设置 `CoordinateFactor=10` 后直接显示为 `um`。显示缩放与物理
单位参数分别保存，不能混用。

- `bin1d` 不含 y 范围，必须显式设置 `TransverseWidth`。
- `bin2d` 默认对完整 y 宽度统计；设置 `SliceCentersY` 和
  `SliceWidthsY` 可同时得到多个 y 位置/宽度的 mass-x 曲线。
- `SphDimension=3` 时单粒子质量改为 `rho0*d0^3`，并必须用
  `OutOfPlaneWidth` 指定 z 方向统计宽度；最终除以 `Wy*Wz`。
- 任意分片边界穿过 y 网格时，程序按网格与分片的几何重叠比例分配该格子的
  `Ncount`，并用实际覆盖宽度归一化。
- 局部面密度表示每个 x chunk 的质量贡献，累计值直接沿 x 求和，不再除以
  x chunk 宽度。
- 该模型假定所选粒子具有相同的初始质量。多材料或自适应变质量结果若只有
  总 `Ncount`，无法从 chunk 文件恢复各材料的精确质量，需分别输出各材料计数。

SPH 文件中已经输出的 `c_rho`、`c_temp`、`c_damagedPress` 等原始场直接读取，
不使用 `dV`。`dV` 参数只属于 MD 守恒量的派生分析路径；压力、密度和应力等
按公式需要正的 `dV`，速度/温度等则由对应 MD 守恒量公式计算。无论哪种情况，
SPH 原始场和 SPH mass-x 都不会进入 MD 的 `dV` 路径。

## Current SPID chunk compatibility

The reader accepts both legacy LAMMPS-style chunk files and the current SPID
format. Current SPID headers may contain:

```text
# Chunk-averaged data for SPID
# Units microscale
# Name bin1d Kind spatial
# Timestep Number-of-chunks Total-count
# Chunk Coord1 Ncount rho
# Time 0.1
```

- Header length is not fixed; `Units`, `Name`, and `Kind` metadata are kept.
- Embedded `# Time` values support `SelectBy=Time` without a Slurm file.
- `spatial` output supports automatic 1D, 2D, and 3D coordinate detection.
- `field` mass-v output automatically uses `Coord1`; legacy files keep
  `Chunk` as the preferred automatic coordinate.
- `cluster` output accepts both SPID `x/y/z` and legacy `c_x/c_y/c_z` names.
- When `UseFileUnitMetadata=true`, known SPID unit systems are converted to
  the displayed mass-v and mass-x units. Metadata-free files keep the existing
  manual factors.

## 图形界面

- “Calculation parameters”随分析类型动态切换，并为枚举值提供下拉框、逻辑值提供复选框。
- “Visualization”默认显示“全部视图”仪表板，也可从下拉框切换到某个单独视图。
- “绘图数据”页可按视图查看图中实际使用的 X/Y/Z 或分布数据；大表按每页 500
  行显示，可复制所选行、复制完整数据或将当前视图直接导出为 CSV。
- 默认勾选 `Separate figures after run`：每个可用视图还会打开为独立 MATLAB 图窗；取消勾选可只保留嵌入式仪表板，也可随时点击 `Open all figures`。
- 主界面提供中文工作流提示和“示例数据”菜单，可一键载入一维、二维、颗粒度、SPH mass-x 一维/二维分片、mass-v 或孔隙网络案例，并自动设置匹配的分析类型与关键参数。
- 快捷键：`F5` 运行、`Esc` 请求取消、`Ctrl+O` 浏览数据、`Ctrl+E` 加载当前分析类型的示例。
- 不同分析拥有各自视图，例如场图、直方图、X/Y 均值剖面、微分/累积分布、组分标签、孔径分布和孔隙率方向剖面。
- `UpdateMode=replace` 会覆盖已有图层；`UpdateMode=overlay` 会在对应视图中叠加不同算例或时刻。
- “Plot conditions”可修改坐标范围、色限、配色、线宽、点大小、网格、色条、等比例轴和对数轴。
- “Output conditions”支持 MAT、摘要 CSV、分析明细 CSV、逐视图绘图数据 CSV、
  PNG、FIG、PDF、运行清单和自动导出。设置 `SavePlotDataCSV=true` 会为每幅
  可用结果图输出一份与绘图完全一致的数值表。
- 配置文件会按参数名升级，新增选项不会导致旧配置错位。
- “工具”菜单提供“检查当前输入与参数”“系统诊断”和“打开运行日志”，可在正式计算前发现文件列缺失、表头冲突和参数错误。

## 稳定性与可诊断性

- 界面和脚本请求共用同一套参数目录与语义校验，未知参数、重复参数、非法范围会在读取大文件前终止。
- 长文件建立索引和回退逐行读取期间支持取消；损坏或过期的 `*.stepidx.mat` 会自动重建，缓存采用临时文件原子替换。
- 输入预检验证首个数据块、数据列数及 MATLAB 字段名冲突，避免两个原始列名归一化后相互覆盖。
- 导出先在隔离临时目录中完整生成，再统一提交到目标目录；任一格式失败时不保留半成品。
- 默认输出目录为项目下的 `outputs/`，不会再把时间戳文件散落到项目根目录；该目录已加入 `.gitignore`。
- 运行日志位于 MATLAB 用户配置目录 `POST_DATA2/logs/postdata.log`，超过 5 MB 自动轮换为 `.1`。

脚本中可直接检查安装环境：

```matlab
postdata_startup
report = pd_system_diagnostics();
disp(report)
```

## 项目结构

```text
POST_DATA2/
  postdata_app.m              GUI 启动入口
  postdata_run.m              唯一分析调度入口
  postdata_startup.m          MATLAB 路径初始化
  src/
    app/                      界面和状态编排
    core/                     请求、参数目录、验证、配置升级
    io/                       chunk/Slurm 读取、索引和预检
    analysis/                 数值分析
      common/                 共享统计与分布算法
      network2d/              二维网络专用模块
    plot/                     多视图渲染
    export/                   数据和图形导出
  fixtures/sample/            最小回归数据
  fixtures/generated/         一维、二维和团簇较大测试数据
  examples/quick_start.m      脚本示例
  examples/generate_large_test_data.m  较大数据生成器
  examples/run_large_verification.m   批量运行与图形验证
  tests/test_postdata.m        自动化测试
```

更详细的模块边界见 [ARCHITECTURE.md](ARCHITECTURE.md)。

## 测试

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA2')
postdata_startup
results = runtests('tests');
assert(numel(results) > 0, 'No tests were discovered.');
assertSuccess(results);
```

测试会生成读取索引 `*.stepidx.mat` 以加速后续访问；该类缓存已加入 `.gitignore`，发布或交付前可直接删除。

重新生成并运行较大测试案例：

```matlab
cd('C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA2')
addpath('examples')
paths = generate_large_test_data();
report = run_large_verification();
```

默认验证图和摘要写入 `outputs/verification/`，该目录是可随时
重建的运行输出，不纳入源码和分发包。

脚本中也可以直接取得某幅图的完整数据：

```matlab
viewData = pd_result_plot_data(result, 'cumulative');
disp(viewData.ColumnNames)
disp(viewData.Values)
pd_write_plot_data_csv('mass_x_cumulative.csv', viewData);
```
