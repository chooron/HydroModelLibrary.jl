3.1 降水分配 (Precipitation Partitioning)

一般容量控制公式


$$R_{int} = \theta_{rain} \cdot R$$

$$S_{int} = \theta_{snow} \cdot S$$

$$LAI = (1-s) \cdot LAI_{max} \cdot f_{LAI}(m)$$

$$SAI = (1-s) \cdot \beta \cdot h_{veg}$$

$$h_{veg} = h_{max} \cdot f_{veg}(m)$$

$$C = C_{max} \cdot \frac{LAI}{LAI_{max}}$$

基于LAI的线性方法 (PRECIP_ICEPT_LAI)


$$\theta_{rain} = \alpha_{rain} \cdot (LAI+SAI)$$

$$\theta_{snow} = \alpha_{snow} \cdot (LAI+SAI)$$

基于LAI的指数方法 (PRECIP_ICEPT_EXPLAI)


$$\theta_{rain} = 1 - \exp(-0.5(LAI+SAI))$$

$$\theta_{snow} = 1 - \exp(-0.5(LAI+SAI))$$

3.2 下渗/径流分配 (Infiltration / Runoff partitioning)

土壤存储一般公式


$$\phi_{max} = Hn(1-SF)$$

$$\phi_{tens} = \phi_{max}(S_{fc}-S_{wilt})$$

$$\phi_{fc} = \phi_{max}S_{fc}$$

推理方法 (INF_RATIONAL)


$$M_{inf} = R \cdot (1-P_{c})$$

SCS方法 (INF_SCS)


$$M_{inf} = R \left( 1 - \frac{(R-0.2S)^2}{R+0.8S} \right)$$

$$S = \frac{25400}{CN} - 254$$

显式与简单Green Ampt方法 (INF_GREEN_AMPT / INF_GA_SIMPLE)


$$M_{inf} = \min \left( R, k_{sat} \left( 1+\frac{|\psi_{f}|(\phi_{max}-\phi_{soil})}{F} \right) \right)$$

VIC方法 (INF_VIC)


$$M_{inf} = R \cdot K_{1} \left( \gamma\alpha z_{max}+z_{min}-\frac{\phi_{soil}}{\phi_{max}} \right)^{\gamma}$$

$$K_{1} = ((z_{max}-z_{min})\alpha\gamma)^{-\gamma}$$

VIC/ARNO方法 (INF_VIC_ARNO)


$$M_{inf} = R \cdot \left( 1-\left( 1-\frac{\phi_{soil}}{\phi_{max}} \right)^{b} \right)$$

HBV方法 (INF_HBV)


$$M_{inf} = R \cdot \left( 1-\left( \frac{\phi_{soil}}{\phi_{max}} \right)^{\beta} \right)$$

PRMS方法 (INF_PRMS)


$$M_{inf} = R \cdot \left( 1-F_{sat}^{max}\min \left( \frac{\phi_{soil}}{\phi_{tens}},1 \right) \right)$$

UBC流域模型方法 (INF_UBC)


$$M_{inf} = R \cdot (1-b_{2})$$

$$b_{2} = b_{1}+(1-b_{1}) \cdot FF$$

$$b_{1} = F_{imp} \cdot 10^{\left( -\frac{\phi_{max}-\phi_{soil}}{PoAGFN} \right)}$$

$$FF = 1+ \frac{\log(\phi_{pond}/VoFLAX)}{\log(VoFLAX/1800)}$$

$$M_{perc} = \min(M_{max}^{perc}, R-M_{inf}) \cdot (1-b_{2})$$

$$M_{int} = (R-M_{inf}-M_{perc}) \cdot (1-b_{2})$$

$$M_{run} = b_{2} \cdot R$$

GR4J方法 (INF_GR4J)


$$M_{inf} = \phi_{max} \cdot \left( \frac{\alpha \cdot \left( 1-\left( \frac{\phi_{soil}}{\phi_{max}} \right)^{2} \right)}{1+\alpha\phi_{soil}} \phi_{max} \right)$$


(其中 $\alpha=\tanh(\phi_{pond}/\phi_{max})$)

HMETS方法 (INF_HMETS)


$$M_{inf} = R \cdot \left( 1-\alpha\cdot\frac{\phi_{soil}}{\phi_{soil}^{max}} \right)$$

AWBM方法 (INF_AWBM)


$$M_{excess} = \max\left(a_{1}R-\frac{\max(a_{1}\phi_{1}^{max}-\phi_{1},0)}{\Delta t},0\right) + \max\left(a_{2}R-\frac{\max(a_{2}\phi_{2}^{max}-\phi_{2},0)}{\Delta t},0\right) + \max\left(a_{3}R-\frac{\max(a_{3}\phi_{3}^{max}-\phi_{3},0)}{\Delta t},0\right)$$

$$M_{inf} = R - M_{excess}$$

$$M_{runoff} = (1-BFI) \cdot M_{excess}$$

$$M_{toGW} = BFI \cdot M_{excess}$$

3.3 基流 (Baseflow)

恒定基流 (BASE_CONSTANT)


$$M_{base} = M_{max}$$

线性存储 (BASE_LINEAR)


$$M_{base} = k\phi_{soil}$$

解析线性存储 (BASE_LINEAR_ANALYTIC)


$$M_{base} = \phi_{soil}\cdot(1-\exp(-k\Delta t))/\Delta t$$

非线性存储 (BASE_POWER_LAW)


$$M_{base} = k\phi_{soil}^n$$

VIC基流方法 (BASE_VIC)


$$M_{base} = M_{max}\left( \frac{\phi_{soil}}{\phi_{max}} \right)^{n}$$

GR4J基流方法 (BASE_GR4J)


$$M_{base} = \frac{\phi_{soil}}{\Delta t}\cdot\left( 1-\left( 1+\left( \frac{\phi_{soil}}{\phi_{ref}} \right)^{4} \right)^{\frac{1}{4}} \right)$$

Topmodel基流方法 (BASE_TOPMODEL)


$$M_{base} = M_{max}\cdot\frac{\phi_{max}}{n}\cdot\frac{1}{\lambda^{n}}\cdot\left( \frac{\phi_{soil}}{\phi_{max}} \right)^{n}$$

基于阈值的幂律基流方法 (BASE_THRESH_POWER)


$$M_{base} = M_{max}\cdot\left( \frac{\frac{\phi_{soil}}{\phi_{max}}-S_{th}}{1-S_{th}} \right)^{n}$$

基于阈值的存储基流方法 (BASE_THRESH_STOR)


$$M_{base} = K_2 \max(\phi_{soil} - \phi_{th}, 0)$$

3.4 深层渗漏 (Percolation)

恒定渗漏 (PERC_CONSTANT)


$$M_{perc} = M_{max}$$

校正线性渗漏 (PERC_GAWSER)


$$M_{perc} = M_{max}\left( \frac{\phi_{soil}-\phi_{fc}}{\phi_{max}-\phi_{fc}} \right)$$

线性渗漏 (PERC_LINEAR)


$$M_{perc} = k\phi_{soil}$$

幂律渗漏 (PERC_POWER_LAW)


$$M_{perc} = M_{max}\left( \frac{\phi_{soil}}{\phi_{max}} \right)^{n}$$

PRMS渗漏方法 (PERC_PRMS)


$$M_{perc} = M_{max}\left( \frac{\phi_{soil}-\phi_{tens}}{\phi_{max}-\phi_{tens}} \right)^{n}$$

Sacramento渗漏方法 (PERC_SACRAMENTO)


$$M_{perc} = M_{max}^{base} \left( 1+\alpha\left( 1-\frac{\phi_{soil}^{to}}{\phi_{max}^{to}} \right)^{\psi} \right) \left( \frac{\phi_{soil}-\phi_{tens}}{\phi_{max}-\phi_{tens}} \right)$$

GR4J渗漏方法 (PERC_GR4JEXCH / PERC_GR4JEXCH2)


$$M_{pere} = -x_{2} \cdot (\min(\phi_{soil}/x_{3},1.0))^{3.5}$$

3.5 壤中流 (Interflow)

PRMS壤中流方法 (INTERFLOW_PRMS)


$$M_{inter} = M_{max} \left( \frac{\phi_{soil} - \phi_{tens}}{\phi_{max} - \phi_{tens}} \right)$$

3.6 土壤蒸散发 (Soil Evapotranspiration)

未校正蒸发 (SOILEVAP_ALL)


$$M_{evap} = PET$$

线性蒸发-饱和度 (SOILEVAP_HBV / SOILEVAP_TOPMODEL)


$$M_{evap} = PET\cdot \min\left( \frac{\phi_{soil}}{\phi_{tens}},1 \right)$$

VIC土壤蒸发 (SOILEVAP_VIC)


$$M_{evap} = PET\cdot\left( 1-\left( 1-\frac{\phi_{soil}}{\phi_{max}} \right)^{\gamma} \right)$$

线性蒸发-存储 (SOILEVAP_LINEAR)


$$M_{evap} = \min(\alpha\cdot\phi_{soil}, PET)$$

根系分布双层蒸发 (SOILEVAP_ROOT)


$$M_{evap}^{U} = PET\cdot\xi_{U}\cdot \min\left( \frac{\phi_{soil}^{U}}{\phi_{tens}^{U}},1 \right)$$

$$M_{evap}^{L} = PET\cdot\xi_{L}\cdot \min\left( \frac{\phi_{soil}^{L}}{\phi_{tens}^{L}},1 \right)$$

顺序双层蒸发 (SOILEVAP_SEQUEN)


$$M_{evap}^{U} = PET\cdot \min\left( \frac{\phi_{soil}^{U}}{\phi_{tens}^{U}},1 \right)$$

$$M_{evap}^{L} = (PET-M_{evap}^{U})\cdot \min\left( \frac{\phi_{soil}^{L}}{\phi_{tens}^{L}},1 \right)$$

UBCWM方法 (SOILEVAP_UBC)


$$M_{evap} = PET\cdot(1-\beta_{fast})10^{\left( -\frac{\phi_{max}-\phi_{soil}}{\gamma_{e}} \right)}$$

$$\beta_{fast} = F_{imp}\cdot10^{\left( -\frac{\phi_{max}-\phi_{soil}}{\gamma_{a}} \right)}$$

GR4J土壤蒸发 (SOILEVAP_GR4J)


$$M_{evap} = \alpha\phi_{soil}\frac{2.0-\frac{\phi_{soil}}{\phi_{max}}}{1.0+\alpha\left( 1.0-\frac{\phi_{soil}}{\phi_{max}} \right)}$$

安大略作物热量单位方法 (SOILEVAP_CHU)


$$M_{evap} = \frac{CHU}{CHU_{mat}}\cdot PET$$

HYPR土壤/湿地蒸发 (SOILEVAP_HYPR)


$$M_{evap}^{*} = PET\cdot \min\left( \frac{\phi_{soil}}{\phi_{tens}},1 \right)$$

$$F_{p} = F_{max}\cdot\left( \frac{\phi_{dep}}{\phi_{dmax}} \right)^{n}$$

$$M_{evap}^{d} = (1-F_{p})\cdot M_{evap}^{*}$$

$$M_{evap}^{s} = F_{p}\cdot PET_{OW}$$

AWBM蒸发方法 (SOILEVAP_AWBM)


$$M_{evap}^{1} = a_{1}PET$$

$$M_{evap}^{2} = a_{2}PET$$

$$M_{evap}^{3} = a_{3}PET$$

HYMOD2蒸发方法 (SOILEVAP_HYMOD2)


$$M_{evap} = K\cdot PET$$

$$K = K_{max}\cdot\left( G+(1-G)\cdot\left( \frac{c^{*}}{c_{max}} \right)^{c} \right)$$

$$c^{*} = c_{max} \cdot \left( 1.0-\left( 1.0-\frac{\phi}{\phi_{max}} \right)^{\frac{1}{b+1.0}} \right)$$

3.7 毛细上升 (Capillary Rise)

HBV模型毛细上升 (CRISE_HBV)


$$M_{crise} = M_{max}^{cr} \left( 1 - \frac{\phi_{soil}}{\phi_{max}} \right)$$

3.9 冠层蒸发 (Canopy Evaporation)

最大冠层蒸发 (CANEVP_MAXIMUM)


$$M_{evap} = PET\cdot F_{c}\cdot(1-f_{s})$$

Rutter冠层蒸发 (CANEVP_RUTTER)


$$M_{evap} = PET\cdot F_{c}\cdot(1-F_{t})\left( \frac{\phi_{can}}{\phi_{cap}} \right)$$

3.10 冠层滴落 (Canopy Drip)

缓排冠层滴落 (CANDRIP_SLOWDRAIN)


$$M_{drip} = \alpha \frac{\phi_{can}}{\phi_{cap}}$$

3.11 洼地蓄水/抽象 (Abstraction)

SCS方法 (ABST_SCS)


$$M_{abst} = \frac{1}{\Delta t}\max\left( f_{SCS}\cdot 25.4 \left( \frac{1000}{CN}-10 \right),\phi_{pond} \right)$$

百分比方法 (ABST_PERCENTAGE)


$$M_{abst} = \alpha M_{pond}$$

PDMROF方法 (ABST_PDMROF)


$$c_{max} = \phi_{max}\cdot(b+1)$$

$$c^{*} = c_{max}\left( 1-\left( 1-\frac{\phi_{dep}}{\phi_{max}} \right)^{\frac{1}{b+1}} \right)$$

$$M_{abst} = \frac{\phi_{max}}{\Delta t} \left[ \left( 1-\frac{c^{*}}{c_{max}} \right)^{b+1} - \left( 1-\frac{c^{*}+\phi_{pond}}{c_{max}} \right)^{b+1} \right]$$

UWFS方法 (ABST_UWFS)


$$D_{min}^{new} = D_{min} - \min(\beta_{ave}P, D_{min})$$

3.12 洼地/湿地蓄水溢流 (Depression/Wetland Storage Overflow)

幂律阈值 (DFLOW_THRESHPOW)


$$M_{dflow} = M_{max}\cdot\left( \frac{\phi_{dep}-\phi_{th}}{\phi_{max}-\phi_{th}} \right)^{n}$$

线性洼地溢流 (DFLOW_LINEAR)


$$M_{dflow} = k_{d}\cdot(\phi_{dep}-\phi_{th})$$

堰式洼地溢流 (DFLOW_WEIR)


$$M_{dflow} = 0.666 \cdot c \cdot \sqrt{2g} \cdot (\max(\phi_{dep}-\phi_{th},0.0))^{1.5}$$

$$c = r_{d}\frac{\sqrt{A}}{A}$$

3.13 洼地/湿地渗漏 (Seepage from Depressions/Wetlands)

线性渗漏 (SEEP_LINEAR)


$$M_{dflow} = k_{seep}\cdot\phi_{dep}$$

3.14 湖泊释放 (Lake Release)

线性释放 (LAKEREL_LINEAR)


$$M_{lrel} = k_{lrel}\cdot\phi_{lake}$$

3.15 开阔水域蒸发 (Open Water Evaporation)

基础方法 (OPEN_WATER_EVAP)


$$M_{evap} = C\cdot PET_{ow}$$

河岸方法 (OPEN_WATER_RIPARIAN)


$$M_{evap} = C\cdot f_{s}\cdot PET_{ow}$$

升尺度湿地充溢方法 (OPEN_WATER_UWFS)


$$M_{evap} = C\cdot f_{d}\cdot PET_{ow}$$

$$D_{new} = D_{old} + M_{evap} \Delta t$$

3.17 积雪平衡 (Snow Balance)

最大液态水容量


$$\phi_{max}^{sl} = SWE\cdot SWI$$

简单融化 (SNOBAL_SIMPLE_MELT)


$$M_{melt} = M_{melt}^{\prime}$$

HBV积雪平衡 (SNOBAL_HBV)


$$M_{refreeze} = K_{a}\cdot \max(T_{f}-T,0)$$

Cema Neige积雪平衡 (SNOBAL_CEMA_NEIGE)


$$M_{melt} = \left( 0.1+0.9\cdot \min\left( \frac{\phi_{SWE}}{S_{Ann}},1 \right) \right)\cdot M^{\prime}$$

HMETS积雪平衡 (SNOBAL_HMETS)


$$M_{rf} = K_{f}\cdot(T_{rf}-T_{di})^{f}$$

$$SWI = \max(SWI_{min}, SWI_{max}(1 - aM_{cumul}))$$

3.18 积雪升华 (Snow Sublimation)

Kuzmin方法 (SUBLIM_KUZMIN)


$$M_{subl} = (0.18+0.098\cdot v_{ave})\cdot(P_{sat}-P_{vap})$$

中央内华达方法 (SUBLIM_CENTRAL_SIERRA)


$$M_{subl} = 0.0063\cdot(h_{w}\cdot h_{v})^{-\frac{1}{6}}\cdot(P_{sat}-P_{vap})\cdot v_{ave}$$

体空气动力学方法 (SUBLIM_BULK_AERO)


$$Q_{e} = \rho_{a}\lambda_{s}\cdot C_{E}\cdot v_{ave}\cdot\frac{0.622(P_{sat}^{snow}-P_{ave})}{P_{a}}$$

$$C_{E} = \frac{\kappa^{2}}{\ln(z_{ref}/z_{0})\ln(z_{ref}/z_{0}^{e})}$$

$$M_{subl} = \frac{Q_{E}}{\lambda_{s}\rho_{w}}$$

3.19 积雪重冻结 (Snow Refreeze)

度日法 (FREEZE_DEGREE_DAY)


$$M_{frz} = K_f \max(T_f - T_a, 0)$$

3.20 积雪反照率演变 (Snow Albedo Evolution)

UBC流域模型方法 (SNOALB_UBC)


$$M_{snalb} = -\alpha\cdot\frac{1-K}{\Delta t}+\frac{(\alpha_{max}-\alpha)}{\Delta t}\min\left( \frac{SN}{SN_{alb}},1 \right)$$

CRHM Essery方法 (SNOALB_CRHM_ESSERY)


$$M_{snalb} = -\beta \quad (\text{if } T_{snow} < 0)$$

$$M_{snalb} = -\beta_{2}\cdot(\alpha-\alpha_{min}) \quad (\text{if } \alpha < \alpha_{b})$$

$$M_{snalb} = (\alpha_{max}-\alpha)\cdot \min(S/S_{thresh},1.0)\cdot\Delta t$$

Baker方法 (SNOALB_BAKER)


$$\alpha = 0.9-0.0473\cdot A_{snow}^{0.1}$$

3.22 冰川释放 (Glacier Release)

线性存储 (GRELEASE_LINEAR_STORAGE)


$$M_{grelease} = -K\phi_{glac}$$

解析线性存储 (GRELEASE_LINEAR_ANALYTIC)


$$M_{grelease} = \frac{\phi_{glac}}{\Delta t}(1-\exp(-K\Delta t))$$

HBV-EC方法 (GRELEASE_HBV_EC)


$$M_{grelease} = -K^{*}\phi_{glac}$$

$$K^{*} = K_{min}+(K-K_{min})\exp(-AG(SN+SN_{liq}))$$

3.24 湖泊结冰 (Lake Freezing)

基础方法 (LFREEZE_BASIC)


$$\frac{dT}{dt} = -\left( 1-\min\left( \frac{\phi_{SWE}}{B},1 \right) \right)M_{melt}$$

3.25 作物热量单位演变 (Crop Heat Unit Evolution)

安大略方法 (CHU_ONTARIO)


$$CHU = 3.33 (T_{max}-10) - 0.084 (T_{max}-10)^2$$

$$CHU_{n} = 1.8\cdot(T_{min}-4.4)$$

$$CHU_{new} = CHU_{old} + 0.5 (CHU_n + CHU)$$

3.26 卷积 (Convolution)

通用卷积积分


$$M_{conv} = \int_{0}^{\infty}UH(\tau)I(t-\tau)d\tau$$

GR4J传递函数 1 (CONVOL_GR4J_1)


$$UH(t) = \begin{cases} \frac{5}{2x_{4}}\left( \frac{t}{x_{4}} \right)^{\frac{3}{2}} & t \le x_{4} \\ 0 & t > x_{4} \end{cases}$$

GR4J传递函数 2 (CONVOL_GR4J_2)


$$UH(t) = \begin{cases} \frac{5}{kx_{4}}\left( \frac{t}{k_{4}} \right)^{\frac{3}{2}} & t \le x_{4} \\ \frac{5}{x_{4}}\left( 2-\frac{t}{x_{4}} \right)^{\frac{3}{2}} & x_{4} < t \le 2x_{4} \\ 0 & t > 2x_{4} \end{cases}$$

Gamma传递函数 1 & 2 (CONVOL_GAMMA, CONVOL_GAMMA2)


$$UH(t) = \frac{1}{t}\frac{(\beta t)^{a}}{\Gamma(a)}\exp(-\beta t)$$