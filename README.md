# POST_DATA2

面向 LAMMPS chunk、cluster、mass-v、mass-x 和二维孔隙网络数据的 MATLAB 后处理软件。项目仅保留一套模块化实现，兼容 MATLAB R2016b，不依赖 App Designer，也不使用 `+package` 目录。

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
request.makePlots = false;
result = postdata_run(request);
```

## 分析类型

- `chunk`：一维或二维原始场/派生场。`ChunkDim=auto` 会根据 `Coord2/c_y` 是否存在自动识别维度，避免一维文件误报缺少 `Coord2`。
- `cluster`：团簇尺寸、概率分布、CDF、拟合和位置分箱统计。
- `vx`：mass-v 面密度分布，默认从最大速度向最小速度累积。
- `massx`：mass-x 空间分布，固定从最大坐标向最小坐标累积；累计曲线按横坐标从小到大显示时保证单调递减，纵轴单位为 `mg/cm^2`。负面密度默认裁剪为零，也可设为报错。
- `network2d`：二维孔隙/基体几何、拓扑、连通性、方向剖面和演化统计。

`massx` 可设置坐标列、坐标换算系数、坐标范围、坐标名称/单位、累积方向和要统计的 `ArealDensity` 列。若密度列留空，程序自动发现所有以 `ArealDensity` 结尾的列；存在 `mass1ArealDensity` 与 `mass2ArealDensity` 时，`massArealDensity` 会按二者之和重建。

## 图形界面

- “Calculation parameters”随分析类型动态切换，并为枚举值提供下拉框、逻辑值提供复选框。
- “Visualization”默认显示“全部视图”仪表板，也可从下拉框切换到某个单独视图。
- 默认勾选 `Separate figures after run`：每个可用视图还会打开为独立 MATLAB 图窗；取消勾选可只保留嵌入式仪表板，也可随时点击 `Open all figures`。
- 主界面提供中文工作流提示和“示例数据”菜单，可一键载入一维、二维、颗粒度、mass-x、mass-v 或孔隙网络案例，并自动设置匹配的分析类型与关键参数。
- 快捷键：`F5` 运行、`Esc` 请求取消、`Ctrl+O` 浏览数据、`Ctrl+E` 加载当前分析类型的示例。
- 不同分析拥有各自视图，例如场图、直方图、X/Y 均值剖面、微分/累积分布、组分标签、孔径分布和孔隙率方向剖面。
- `UpdateMode=replace` 会覆盖已有图层；`UpdateMode=overlay` 会在对应视图中叠加不同算例或时刻。
- “Plot conditions”可修改坐标范围、色限、配色、线宽、点大小、网格、色条、等比例轴和对数轴。
- “Output conditions”支持 MAT、摘要 CSV、明细 CSV、PNG、FIG、PDF、运行清单和自动导出。
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

结果图和数值摘要输出到 `verification/`。
