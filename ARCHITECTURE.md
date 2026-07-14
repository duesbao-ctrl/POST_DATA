# POST_DATA non-GUI architecture

依赖方向固定为：

```text
run_analysis -> postdata_run -> core/io/analysis
postdata_run -> plot (仅在请求生成独立图时统一样式)
plot -> core
export -> core
analysis -> io
```

- `postdata_run` 是唯一调度器；`run_analysis` 仅转换旧参数到版本化请求。
- `src/core` 是参数、验证、选择方式和结果契约的唯一事实来源；脚本和界面请求均
  必须通过同一参数目录校验。
- `src/analysis` 只负责数值，`network2d` 的网格、几何、拓扑、连通性、方向剖面
  和演化均拆分为独立模块。
- `src/plot` 将结果渲染到调用方坐标轴，并统一应用可编辑论文样式。
- `src/io` 集中处理表头歧义、预检、可取消索引读取和原子缓存，分析模块不得直接
  依赖侧车缓存文件。
- `src/export` 独立处理数据表、位图、矢量 PDF/FIG 和可复现清单，并通过临时区
  完整生成后统一提交。
- 本项目禁止引入 GUI 状态或复制分析公式；界面编排只存在于 `POST_DATA2/src/app`。

增加功能时先扩展参数目录与结果契约，再实现分析，最后增加渲染/导出与测试。
所有生产代码必须保持 MATLAB R2016b 兼容。
