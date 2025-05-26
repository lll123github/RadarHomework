%RCS建模 Swerling2模型
tgtmodel = 'Swerling2';
f=2e9; % 工作频率 S波段
antenna = phased.IsotropicAntennaElement('BackBaffled',true);
radiator = phased.Radiator('OperatingFrequency',f,'Sensor',antenna);
RCS=6;
target = phased.RadarTarget('MeanRCS',RCS,'OperatingFrequency',f,...
    'Model',tgtmodel);

sigma = target.MeanRCS; % 目标雷达散射截面

p_t=40*10e3; % 发射功率

lambda=1/f;% 波长
theta_B=1.1/360*pi;
phi_B=5.7/360*pi; % 波束宽度
G_D=4*pi/(theta_B*phi_B);
rau=0.7; % 天线效率
G=G_D*rau; % 天线增益

SN_min_db=13; % 信噪比
SN_min=10^(SN_min_db/10); % 最小信噪比
L_s=10^(L_s_dB/10); % 系统损耗
% L_a=10^(L_a_dB/10); % 大气损耗
L_a=1; % 大气损耗

p_r=(p_t*G^2*lambda^2*sigma)/((4*pi)^3*L_s*L_a); % 接收功率

PRF=1e3; % 脉冲重复频率
tau=1*10e-6; % 脉冲宽度
% f_s=1e11; % 采样频率
% B_n=f_s; % 噪声带宽
T=6/60; % 脉冲重复周期 每分钟六转
T_0=290; % 参考温度
F_n_dB=3.5; % 噪声系数(dB)
F_n=10^(F_n_dB/10); % 噪声系数
k=1.38e-23; % 玻尔兹曼常数

R_max=(p_r*tau*PRF*T/(F_n*k*T_0*SN_min))^(1/4); % 最大探测距离