# Skim → fake rate → Z+X

独立的 Python/PyROOT 分析程序，读取 HHbbZZ 的 `Events` skim tree。
支持每个 era、每个轻子通道及所有 era 合并结果；所有输入文件、样本分类、
归一化、额外选择、branch 表达式和分 bin 均由 YAML 定义。

## 必须先补齐 pass probe

旧版 skim 只有 `passedZXCR2P1FSelection`（fail probe），不能单独确定
`f = N(tight) / N(loose)`。本次同步给上游加入
`passedZXCR3P0FSelection`：同一个 Z1+l 选择下，第三个轻子通过 tight。
2P1F 与 3P0F 都要求恰好三个 loose 轻子、tight Z1、MET < 25 GeV、
|mZ1−mZ| < 7 GeV，以及至少两个通过选择及候选轻子重叠清理的 jet。
3P0F 和 2P1F 共用最接近 Z 质量的 Z1 选择。

**从数据测量 fake rate，需要用更新后的上游重新产生 data 和 MC skim。**
旧文件仍可用于原有三个 CR 的绘图、或读取外部 fake-rate map 后做 OS 估计；
不能从已经丢弃的事件恢复 3P0F。程序遇到缺失的 branch 会报错。

输入约定：

| 用途 | Branch |
|---|---|
| FR fail | `passedZXCR2P1FSelection` |
| FR pass | `passedZXCR3P0FSelection` |
| 两个 fail 的应用区 | `passedZXCR2P2FSelection` |
| 一个 fail 的应用区 | `passedZXCR3P1FSelection` |
| probe 运动学 | `pTL3`, `etaL3`, `pTL4`, `etaL4` |
| probe 类型及状态 | `ZXCRPdgIdL3/4`, `ZXCRTightL3/4` |
| jet 数 | `njets_pt30_eta4p7`（名称沿用上游，实际阈值取决于 skim 配置） |

三轻子样本只读 L3；四轻子样本根据 tight 标志找 fail，不依赖其排序。
已有四轻子 SR 事件不会进入估计。当前 skim 没有 SS 控制区，故实现的是
**OS 方法**，不声称复现参考分析的 OS/SS 联合估计。

## 配置与运行

在带 PyROOT 的环境（例如相应 CMSSW 环境）使用 Python >= 3.9，并安装 PyYAML。
无需 uproot、pandas 或 matplotlib。所有图使用 ROOT 生成。

```bash
cp ZXEstimation/config.example.yaml ZXEstimation/config.yaml
export SKIM_DIR=/absolute/path/to/skimmed/files
# 编辑 config.yaml：填写 lumi_fb、文件路径、xsec_pb 和 sum_gen_weights
python ZXEstimation/run.py --config ZXEstimation/config.yaml --stage all
```

PowerShell 可用 `$env:SKIM_DIR = 'E:\path\to\skims'`。
配置中的 null 是刻意保留的必填项，不提供虚构的截面、亮度或生成事件数。
相对路径相对于 YAML 目录；也可使用显式 `root://.../file.root` URL。
本地文件支持 glob；远程 URL 不支持 glob。复制样本项即可扩展 era 和数据集。

```bash
python ZXEstimation/run.py --config ZXEstimation/config.yaml --stage plots
python ZXEstimation/run.py --config ZXEstimation/config.yaml --stage fake-rate
python ZXEstimation/run.py --config ZXEstimation/config.yaml --stage estimate
```

`estimate` 读取此前生成的 `output/fake_rates.root`；不会重新求解。
`all` 先画 CR 图，再求 FR，最后估计。若 FR 稀疏导致应用失败，CR 图仍然可用。
用旧 skim 时，设置 `fake_rate.source: external`，配置外部文件，
并把 `plotting.regions` 中的 `3P0F` 删除；不要把 3P0F branch 替换为其他事件类别。
外部文件需提供每个 era/flavor 的 TH2 `rate`（bin error 为 FR 误差）和 TH2 `valid`
（有效 bin = 1），路径可用 `{era}`、`{flavor}` 模板，坐标分 bin 必须与 YAML 一致。

额外选择用 ROOT `TTreeFormula` 标量表达式配置：

- `selection.cut`：所有阶段共用；`min_jets` 不允许小于 2。
- `samples[].selection`：数据流 trigger routing 等样本特定选择。
- `regions.<name>.cut`：该区域的额外选择。例如可给两个应用区同时加相同 m4l 窗口。
- `channels`：任意命名的通道选择。默认包含 inclusive、4e、4mu、2e2mu、2mu2e。
- `variables`：表达式、边界、label、适用 CR；`estimate: true` 同时生成 Z+X 模板。

**输入 skim 的 lepton ID、FSR、jet、trigger 和 era 定义必须与 FR 测量一致。**
示例图的默认变量含 Z/4l 质量和 pT、jj 质量、轻子及 jet 的 pT/η、jet 数；
其他标量变量可直接在 YAML 加条目。jet 1/2 沿用上游 b-tag 排序。

## 权重与扣除

MC 的默认归一化为：

`w = [1000 × lumi_fb × xsec_pb × k_factor × filter_efficiency / sum_gen_weights] × weight_expression`

`sum_gen_weights` 必须是与这些文件对应的**完整生成样本、skim 前**的带符号权重和；
不能用 skim 通过事件数或 skim 的 Weight 总和替代。文件分片和扩展样本需使用匹配的
联合分母；不要把不同截面的 qqZZ、ggZZ 等数据集放进同一个 normalization 项。
它们可各有一个 sample，但共用 `process: ZZ`。

如果 `weight_expression` 已含亮度归一化，改用：

```yaml
normalization: {mode: preweighted, scale: 1.0}
```

数据权重固定为 1，并默认跨数据文件按 `(era, run, luminosityBlock, event)` 去重。
去重不替代 primary-dataset trigger routing；需要时用各 sample 的 selection 配置。

FR 对每个 era、|PDG ID|、pT 和 |η| bin 分别求解：

`T = data(3P0F) − prompt_MC(3P0F)`

`F = data(2P1F) − prompt_MC(2P1F)`

`f = T/(T+F)`

示例 FR 扣除 WZ、ZZ，列表可改。DY、TTbar 的非 prompt 成分不作为 FR 污染扣除。
T/F 使用独立样本的 sumw2，误差为
`Var(f) = [F² Var(T) + T² Var(F)]/(T+F)^4`，保留分子分母相关性。
测量中 T 或 F 非正的 bin 标为无效。应用时使用无效 bin 会失败，需合并分 bin 或
提供经过验证的外部测量，不会自动设 f=0。pT 上溢默认合并到最后一个 bin；
pT 下溢、η 越界始终报错。无效 bin 在 FR 图中灰色显示。

OS 外推使用 `r=f/(1−f)`：

`N(ZX) = Σ[data,3P1F] w r − Σ[ZZ,3P1F] w r − Σ[data,2P2F] w r1 r2`

2P2F 是减项，避免在 3P1F 外推中重复计算两假轻子过程。
WZ、TTbar、DY 在四轻子应用区属于待估计的可约背景，默认不作为 ZZ 污染一起扣掉。
这与 [CMS HIG-16-041，7.2 节，式 (4)](https://arxiv.org/html/1706.09936#S7.SS2)
的 OS 组合一致。参考本地 `yhbbzz/zxcr.py` 的选择和 per-era 流程；本地副本缺少
`studies/zxfr_an2023_157/zxfr_core.py`，本框架无需导入该模块。

`estimation.subtract_single` 与 `subtract_double` 可分别定义污染扣除；默认后者为空。
若启用 double 的 prompt 扣除，它对最终式子的贡献符号为正。

## 输出

默认写入 YAML `output` 目录：

- `control_regions.root`：`era/CR/channel/variable/{data,ZZ,WZ,TTbar,DY,total_mc}`。
- `control_regions.json`：不依赖 histogram 范围的每区域产额、sumw2、去重计数、无定义 ratio bin。
- `plots/era/CR/channel/*.png,pdf`：data 点、四类 MC 的堆叠、MC 统计误差带、data/MC ratio。
- `fake_rates.root`：`era/11`、`era/13` 下的 rate、valid、原始和扣除后的 T/F。
- `fake_rates.json`：逐 bin 的 T/F、FR、误差和有效状态；`plots/fake_rates` 为二维 FR 图。
- `zx_estimation.root`：`era/channel/variable` 和 `combined/channel/variable`。
  `zx_stat` 保存应用区统计误差，`zx_total` 再加入 FR 统计误差；
  `covariance_stat_plus_fr` 保存 bin 间协方差；四个 component 保留各项带符号贡献。
  `variations/fr_<era>_<flavor>_<ptbin>_<etabin>Up/Down` 为逐 FR bin 的一阶误差变化。
- `zx_estimation.json`：各通道总产额、误差、分项和负 bin 列表。`yield` 使用独立的一格直方图，
  因此产额不受所画变量范围影响。不要将 inclusive 和各独立末态再次相加。
- `manifest_<stage>.json`、`config_<stage>.yaml`：本次实际文件列表、归一化及执行状态。

同一 FR bin 被多个事件或两个 probe 使用时按同一 nuisance 传播，era 间独立。
这里的误差为统计与 FR 统计误差的一阶传播；不包含 lumi、截面、ID/SF、fake composition
等额外系统误差，也没有擅自套用参考分析的最终模型参数。
不会截断负 MC 权重或负 Z+X 结果；MC 正、负贡献分别堆叠。ratio 仅对 MC 总量 > 0
的 bin 定义，data 使用 Garwood 68.27% Poisson 误差，灰带为 MC 统计相对误差。
默认将图的上下溢折入边缘；关闭 fold 后 flow 仍写在 ROOT 文件中但不显示在图内。

尚未提供实际 skim 路径和归一化数值，因此代码交付不包含实际 FR 数值或物理产额。
