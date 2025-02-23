** 其他语言版本：[English](README.md)

# rBTR 

rBTR 是 Broadleaved Tree-Ring (BTR) 模型的R包，详细信息可以参看论文：

> Zhao, B. *et al.* A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR). *Tree Physiol.* **44**, tpae127 (2024).[DOI:[10.1093/treephys/tpae127](http://dx.doi.org/10.1093/treephys/tpae127) ]

## 安装R包

可以使用以下代码安装 rBTR 包:

```R
# devtools::install_github("NEFU-TreeRingLab/rBTR")
```
---
## 使用方法

#### 运行模型

模型包含两个主要函数:
1. `btr()`: 单线程运行
2. `btr_parallel()`:支持并行运算
两个函数使用方法一致，在进行短序列模拟时推荐使用`btr()`, 长序列模拟时 `btr_parallel()` 的运行速度更有优势。

```R
# 运行BTR
# 单核版本
# btr(clim = Clims, parameters = BPparam, age = BPage,writeRes = F,Named = 'TestData' )

# 并行运算版本
# btr_parallel(clim = Clims, parameters = BPparam, age = BPage,  writeRes = F,Named = 'TestData' )

# Note:
# clim 为模型输入的气象数据  data.frame
# paramters 为模型参数  data.frame
# age 为模型输入的年龄趋势 data.frame
# writeRes 布尔值，为 TURE 时会模型会在工作目录建立一个文件夹，自动保存模拟结果到'res_(running time)"的文件夹中
# Named 文本，当 "writeRes" 为 TURE 时会模型会在工作目录建立一个以Named中文本为文件名的文件夹并自动保存模拟结果。
```

- 模型参数化方法可以查阅R包 `BTRtools`[Link:BTRtools](https://github.com/NEFU-TreeRingLab/BTRtools)

#### 组织输入数据

- 气象数据

  模型输入的气象数据需要组织为以下格式
  
  | Year | Month | Day  | TEM   | VPD   | soilM |
  | ---- | ----- | ---- | ----- | ----- | ----- |
  | 2020 | 1     | 1    | -23.7 | 0.029 | 0.31  |
  
  模型输入中 ‘***VPD***’ 可以使用相对湿度 ‘***RH***’ 代替
  土壤湿度 ‘***soilM***’ 可以通过降雨 ‘***PRE***’ 来进行估计（降雨估算土壤湿度的模块估算效果没有经过细致的测试）
  
  模型输入数据需要进行一些预计算。需要输入参数 ’***lat***‘ 为样地纬度， ’***rootd***‘ 为样地 根区土层深度.
  
  ```R
  # 模型输入数据
  # data(Climdata) 
  # 初步处理输入数据
  # Clims <- Compute_clim(Climdata, lat = 47.13, rootd = 1000 )
  
  # 使用替代气象数据的模型输入
  # data(Climdata.2)
  # 初步处理输入数据
  # Clims <- Compute_clim(Climdata.2, lat = 47.13, rootd = 1000 )
  ```
  
  初步处理输入的模型输入数据为：
  
  | Year | DOY  | TEM    | VPD   | soilM | gE    | Ls    | dL_i   | rootd |
  | ---- | ---- | ------ | ----- | ----- | ----- | ----- | ------ | ----- |
  | 2020 | 1    | --23.7 | 0.029 | 0.31  | 0.532 | 0.697 | 0.0144 | 1000  |
  
  其中，样地日照时长’**LS**‘，根区土层深度 ’**rootd**‘，可以替换为实测值，以获得更有效的模拟。（示例数据已提供在R包中）

#### 年龄趋势

- 年龄趋势输入：由于树木年龄趋势的形成较为复杂，模型并不直接模拟树木生长的年龄趋势，需要将年龄趋势进行输入，年龄趋势的格式为：

  | Year | age  | Tage      | Lage      |
  | ---- | ---- | --------- | --------- |
  | 2000 | 41   | 0.1718823 | 0.9606754 |

  ***Tage*** 为形成层活性的年龄趋势（使用年轮宽度拟合）

  ***Lage*** 为细胞扩大上限的年龄趋势
  
  
  
  （示例数据已提供在R包中）
  
  ``` R
  # 输入的年龄趋势
  # data(BPage)
  ```

#### 模型参数表

- 模型参数较多，因此模型中汇总为一个 data.frame 进行输入。模型参数表已在R包中提供，推荐导出为csv或xlsx文件进行修改。

  ```R
  # 模型参数表
  # data(BPparam)
  # 将参数表保存至本地
  # write.csv(BPparam,"BPparam.csv" )
  # openxlsx::write.xlsx(BPparam,"BPparam.xlsx" )
  ```

  

