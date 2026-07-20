# POST_DATA non-GUI architecture

## Current SPID input contract

`src/io/pd_read_chunk_preamble` and
`src/io/pd_read_next_chunk_frame_header` are the only chunk header/frame
parsers. They accept legacy headers and current SPID `Units`, `Name`, `Kind`,
and per-frame `Time` metadata without assuming a fixed header length.
`pd_spid_unit_info` owns SPID unit conversions. Spatial 1D/2D/3D, field, and
cluster outputs all enter the existing analysis APIs through this shared I/O
contract.

依赖方向固定为：

```text
run_analysis -> postdata_run -> core/io/analysis
postdata_run -> plot (仅在请求生成独立图时统一样式)
plot -> core
export -> core/plot
analysis -> io
```

- `postdata_run` 是唯一调度器；`run_analysis` 仅转换旧参数到版本化请求。
- `src/core` 是参数、验证、选择方式和结果契约的唯一事实来源；脚本和界面请求均
  必须通过同一参数目录校验。
- `src/core/resources` 保存 UTF-8 中文文本，由 `pd_ui_text` 显式解码；
  可执行 `.m` 文件必须保持纯 ASCII，不依赖 Windows 系统代码页。
- `src/analysis` 只负责数值，`network2d` 的网格、几何、拓扑、连通性、方向剖面
  和演化均拆分为独立模块。
- SPH mass-x 独立于 mass-v 的面密度列算法：它只读取 chunk `Ncount`，
  以 `rho0*d0^2`（二维 SPH）或 `rho0*d0^3`（三维 SPH）构造粒子质量，
  再按 y 宽度及可选 z 宽度归一化。chunk 分箱维度与 SPH 物理维度是独立参数；
  bin2d 任意 y 分片使用分片和网格的覆盖比例。显示缩放不参与物理质量换算。
- SPH 原始场和 SPH mass-x 绕过 MD 派生量路径；`dV` 选项只属于 MD 派生量，
  由各 MD 公式决定是否实际使用。
- `src/plot` 中的 `pd_result_plot_data` 是每幅图实际数值的唯一契约；直方图、CDF、
  分箱均值、方向剖面和二维网格只在这里构造一次。`pd_render_result`、界面数据表、
  剪贴板和 CSV 导出均消费同一契约，不得分别重复计算。
- `src/io` 集中处理表头歧义、预检、可取消索引读取和原子缓存，分析模块不得直接
  依赖侧车缓存文件。
- `src/export` 独立处理分析明细、逐视图数据表、位图、矢量 PDF/FIG 和可复现
  清单，并通过临时区完整生成后统一提交。
- 本项目禁止引入 GUI 状态或复制分析公式；界面编排只存在于 `POST_DATA2/src/app`。

增加功能时先扩展参数目录与结果契约，再实现分析，最后增加渲染/导出与测试。
所有生产代码必须保持 MATLAB R2016b 兼容；新增中文文本必须进入
UTF-8 资源目录，并通过 `pd_ui_text` 访问。
