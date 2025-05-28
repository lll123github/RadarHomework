% 参数设置
tgtmodel = 'Swerling2';
f = 2e9; % S波段频率
c= physconst('LightSpeed'); % 光速
lambda = c/f; % 波长计算

% 天线参数
theta_B_deg = 1.1; % 方位波束宽度（度）
phi_B_deg = 5.7;   % 俯仰波束宽度（度）
theta_B = deg2rad(theta_B_deg);
phi_B = deg2rad(phi_B_deg);
rau = 0.7; % 天线效率
G_D_max = (4*pi)/(theta_B*phi_B); % 方向性
G_max = G_D_max * rau; % 天线增益（线性值）
G_dB_max = 10*log10(G_max); % 天线增益（dB）
arrayLength = [7.3 1.52]; % 阵列长度（米）

%天线建模
%计算阵元间距
d=sqrt(arrayLength(1)*arrayLength(2)*theta_B*phi_B/(4*pi)); % 假设阵列为矩形阵列
%计算d和lamda的比值
d_ratio = d/lambda; % 阵元间距与波长的比值‘
% 打印结果
disp(['阵元间距: ' num2str(d) ' 米']);
disp(['阵元间距与波长的比值: ' num2str(d_ratio)]);

% 计算阵元间距和数量（基于均匀矩形阵列模型）
N_az = round(arrayLength(1)/d); % 方位向阵元数
N_el = round(arrayLength(2)/d); % 俯仰向阵元数

theta = linspace(-pi/2, pi/2, 200);
phi = linspace(-pi/2, pi/2, 200);
[thetaGrid, phiGrid] = meshgrid(theta, phi);

% 计算方向图
G = arrayGain(thetaGrid, phiGrid, N_az, N_el, lambda, d, G_max);

% 显示二维热力图
f1=figure(1);
imagesc(rad2deg(theta), rad2deg(phi), 10*log10(G));
xlabel('Azimuth Angle (deg)');
ylabel('Elevation Angle (deg)');
title('2D Gain Pattern (dB)');
colorbar;
axis xy;

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
    numerator = p_t * G_max^2 *( lambda^2 * RCS * n_p * tau);
    R_max = (numerator/denominator)^(1/4);
    % 显示结果
    disp(['脉冲重复频率: ' num2str(PRF/1e3) ' kHz, 最大探测距离: ' num2str(R_max/1e3) ' km']);
    R_max_list(i) = R_max; % 存储结果
end
% 绘图(对数)
f2=figure(2);
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

% 对不同的PRF绘制不同的最大探测距离热力图

fig_index=3; % 图形索引
for i=1:length(PRF_list) 
    PRF = PRF_list(i);
    n_p = round(PRF * T_dwell); % 积累脉冲数

    % 绘图
    
    G= arrayGain(thetaGrid, phiGrid, N_az, N_el, lambda, d, G_max);
    R=radarEquation(p_t, G, lambda, RCS, n_p, tau,F_n, L_s, L_a, SNR_min);
    figure(fig_index);
    imagesc(rad2deg(theta), rad2deg(phi), R);
    xlabel('Azimuth Angle (deg)');
    ylabel('Elevation Angle (deg)');
    title(['Max Detection Range for PRF = ', num2str(PRF/1e3), ' kHz']);
    colorbar;
    fig_index = fig_index + 1; % 更新图形索引
    axis xy;
end



function G = arrayGain(theta, phi, N_az, N_el, lambda, d, G_max)
        % theta, phi: 角度 (rad)
        % N_az, N_el: 方位与俯仰方向阵元数
        % lambda: 波长
        % d: 阵元间距
        % G_max: 最大增益（线性值）

        k = 2*pi/lambda;

        % 相位项
        psi_az = k * d * sin(theta) .* cos(phi);
        psi_el = k * d * sin(phi);

        % 防止除零
        psi_az(psi_az==0) = 1e-12;
        psi_el(psi_el==0) = 1e-12;

        % 归一化阵列因子
        AF_az = sin(N_az * psi_az / 2) ./ (N_az * sin(psi_az / 2));
        AF_el = sin(N_el * psi_el / 2) ./ (N_el * sin(psi_el / 2));

        AF = AF_az .* AF_el;

        % 增益函数
        G = G_max * abs(AF).^2;
end

function R= radarEquation(p_t, G, lambda, RCS, n_p, tau,F_n, L_s, L_a, SNR_min)
    % 雷达方程计算最大探测距离
    % p_t: 发射功率 (W)
    % G: 天线增益 (线性值)
    % lambda: 波长 (m)
    % RCS: 目标雷达散射截面 (m^2)
    % n_p: 积累脉冲数
    % tau: 脉冲宽度 (s)
    % L_s: 系统损耗 (线性值)
    % L_a: 大气损耗 (线性值)
    % SNR_min: 最小信噪比 (线性值)
    T0 = 290; % 噪声温度（K）
    k = 1.38e-23; % 玻尔兹曼常数
    numerator = p_t .* G.^2 .* (lambda^2 * RCS * n_p * tau);
    denominator = (4*pi)^3 * k * T0 * F_n * L_s * L_a * SNR_min;
    
    R = (numerator ./ denominator).^(1/4); % 最大探测距离
end