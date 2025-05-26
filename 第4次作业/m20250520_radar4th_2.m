% 参数设置
tgtmodel = 'Swerling2';
f = 2e9; % S波段频率
lambda = physconst('LightSpeed')/f; % 波长计算

% 天线参数
theta_B_deg = 1.1; % 方位波束宽度（度）
phi_B_deg = 5.7;   % 俯仰波束宽度（度）
theta_B = deg2rad(theta_B_deg);
phi_B = deg2rad(phi_B_deg);
rau = 0.7; % 天线效率
G_D = (4*pi)/(theta_B*phi_B); % 方向性
G = G_D * rau; % 天线增益（线性值）
G_dB = 10*log10(G); % 天线增益（dB）

% 雷达系统参数
p_t = 40e3; % 发射功率（W）
tau = 1e-6; % 脉冲宽度（s）
PRF = 1e3;  % 脉冲重复频率（Hz）
scan_rate_rpm = 6; % 扫描速率（转/分钟）
scan_rate = 360*scan_rate_rpm/60; % 扫描速率（度/秒）
T_dwell = theta_B_deg/scan_rate; % 驻留时间（秒）
n_p = round(PRF * T_dwell); % 积累脉冲数

% 目标参数
RCS = 6; % 目标RCS（m?）
SNR_min_dB = 13; % 检测所需最小SNR（dB）
SNR_min = 10^(SNR_min_dB/10); % 线性值

% 系统损耗
L_s_dB = 3; % 系统损耗（dB）
L_s = 10^(L_s_dB/10);
L_a_dB = 0; % 大气损耗（dB）
L_a = 10^(L_a_dB/10);

% 噪声参数
F_n_dB = 3.5; % 接收机噪声系数（dB）
F_n = 10^(F_n_dB/10);
T0 = 290; % 噪声温度（K）
k = 1.38e-23; % 玻尔兹曼常数

% 雷达方程计算最大探测距离
numerator = p_t * G^2 * lambda^2 * RCS * n_p * tau;
denominator = (4*pi)^3 * k * T0 * F_n * L_s * L_a * SNR_min;
R_max = (numerator/denominator)^(1/4);

% 显示结果
disp(['最大探测距离: ' num2str(R_max/1e3) ' km']);

% 使用Phased Toolbox验证
waveform = phased.RectangularWaveform('PulseWidth',tau,'PRF',PRF);
transmitter = phased.Transmitter('PeakPower',p_t,'Gain',G);
target = phased.RadarTarget('MeanRCS',RCS,'OperatingFrequency',f,'Model',tgtmodel);

% 计算接收信号功率
[~, range_grid] = rangeaziwaveform(waveform);
tx_sig = transmitter(waveform());
tx_sig = radiator(tx_sig, [0;0]);    % 假设正对目标方向
rx_sig = target(tx_sig);             % 目标反射
rx_sig = receiver(rx_sig);           % 接收机处理

% 计算实际SNR
noise_bw = 1/tau;
SNR_actual = pow2db(bandpower(rx_sig)/(10*log10(k*Ts*noise_bw) + Fn + L_total));
disp(['实际SNR: ' num2str(SNR_actual) ' dB']);


% 方向图可视化（可选）
antenna = phased.CosineAntennaElement('FrequencyRange',[f-0.1e9 f+0.1e9],...
    'CosinePower',[1.5 1.5]);
pattern(antenna,f,'Type','powerdb');