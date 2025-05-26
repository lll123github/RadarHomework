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

scan_rate_rpm = 6; % 扫描速率（转/分钟）
scan_rate = 360*scan_rate_rpm/60; % 扫描速率（度/秒）
T_dwell = theta_B_deg/scan_rate; % 驻留时间（秒）


% 目标参数
RCS = 6; % 目标RCS
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

denominator = (4*pi)^3 * k * T0 * F_n * L_s * L_a * SNR_min;


PRF_list=[1e4,2e4, 4e4, 1e5]; % 脉冲重复频率列表
R_max_list = zeros(size(PRF_list)); % 初始化最大探测距离列表
for i = 1:length(PRF_list) 
    PRF = PRF_list(i);
    n_p = round(PRF * T_dwell); % 积累脉冲数
    % 雷达方程计算最大探测距离
    numerator = p_t * G^2 * lambda^2 * RCS * n_p * tau;
    R_max = (numerator/denominator)^(1/4);
    % 显示结果
    disp(['脉冲重复频率: ' num2str(PRF/1e3) ' kHz, 最大探测距离: ' num2str(R_max/1e3) ' km']);
    R_max_list(i) = R_max; % 存储结果
end
% 绘图(对数)
figure;
plot(PRF_list/1e3, R_max_list/1e3, '-o');
set(gca, 'XScale', 'log'); % 设置x轴为对数刻度
set(gca, 'YScale', 'linear'); % 设置y轴为线性刻度
xlabel('脉冲重复频率 (kHz)');
ylabel('最大探测距离 (km)');




% PRF = 1e5;  % 脉冲重复频率（Hz）需要考虑不同参数的结果
% n_p = round(PRF * T_dwell); % 积累脉冲数
% % 雷达方程计算最大探测距离
% numerator = p_t * G^2 * lambda^2 * RCS * n_p * tau;
% R_max = (numerator/denominator)^(1/4);
% % 显示结果
% disp(['最大探测距离: ' num2str(R_max/1e3) ' km']);




% % 方向图可视化（可选）
% antenna = phased.CosineAntennaElement('FrequencyRange',[f-0.1e9 f+0.1e9],...
%     'CosinePower',[1.5 1.5]);
% pattern(antenna,f,'Type','powerdb');

% % 构建theta和phi网格
% [theta, phi] = meshgrid(linspace(-pi/6, pi/6, 100), linspace(-pi/2, pi/2, 100));

% % 构建方向图模型（例如高斯或cos型）
% gain = cos(theta/(theta_B/2)).^2 .* cos(phi/(phi_B/2)).^2;
% gain(gain < 0) = 0;
% gain_dB = 10*log10(gain);
% gain_dB(gain == 0) = -60;

% % 绘图
% figure;
% surf(rad2deg(theta), rad2deg(phi), gain_dB, 'EdgeColor', 'none');
% xlabel('Azimuth (deg)');
% ylabel('Elevation (deg)');
% zlabel('Gain (dB)');
% title('2D Radiation Pattern');
% colorbar;