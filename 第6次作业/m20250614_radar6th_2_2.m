clc; clear;

% 基础参数
fs = 20e6;          % 采样率 10 MHz
T = 20e-6;          % 脉冲宽度 20 us
B = 5e6;            % 带宽 5 MHz
t = 0:1/fs:T-1/fs;  % 时间轴
N = length(t);      % 点数

% LFM信号
k = B/T;  % 斜率
x = exp(1j*pi*k*t.^2);  % 发射LFM信号（单位幅值）
x=chirp(t,0,T,B); % 线性调频信号
% x=x/sqrt(sum(x.^2));%归一化能量为1

% 匹配滤波器（时域匹配）
h = conj(fliplr(x));  % 匹配滤波器核

%接收信号的信号部分
% y_s=conv(x,h);
% y_s=y_s/sqrt(sum(y_s.^2));%归一化

% SNR范围（单位dB）
SNR_dB = -30:1:0;
SNR_lin = 10.^(SNR_dB/10);

% 蒙特卡洛次数
N_MC = 1e4;

% 门限设定
P_FA_target = 1e-2;  % 虚警概率目标



% 初始化统计值
P_D_sim = zeros(size(SNR_lin));
P_FA_sim = zeros(size(SNR_lin));
P_D_target_arr = zeros(size(SNR_lin));
P_FA_target_arr = zeros(size(SNR_lin));

T_thresh_arr=zeros(1,length(SNR_lin)); % 门限值


for idx = 1:length(SNR_lin)
    snr = SNR_lin(idx);
    detect_count = 0;
    false_alarm_count = 0;
    truth_count=0;
    sqrt_snr = sqrt(snr);
    sigma=sqrt(1/snr);
    % 计算门限
    T_thresh=norminv(1 - P_FA_target, 0, 1) * sigma; % 门限值
    % T_thresh = erfcinv(2* P_FA_target) * sqrt_snr; % 门限值
    T_thresh_arr(idx) = T_thresh; % 存储门限值
    for k = 1:N_MC
        % 随机生成接收信号
        ran_num=rand(1);
        % 以50%的概率生成1
        A=floor(2*ran_num);
        truth_count=truth_count+A;
        %信号部分
        n=sigma*randn(1,N);
        y=conv(A*x+n,h,'same');

        % detect_peak = max(abs(y));
        detect_peak=y(ceil(N/2));
        if ran_num>=0.5 && detect_peak > T_thresh
            detect_count = detect_count + 1; % 检测到目标
        elseif ran_num < 0.5 && detect_peak > T_thresh
            false_alarm_count = false_alarm_count + 1; % 虚警
        end
        
    end
    % 使用一次仿真统计直接计算：
    P_D_sim(idx) = detect_count / truth_count;
    P_FA_sim(idx) = false_alarm_count / (N_MC-truth_count);  
    P_D_target_arr(idx)=1-normcdf((T-1)/sigma); % 目标检测概率 %为什么？
end

f1=figure(1);
semilogy(SNR_dB, P_D_sim, 'ro-', 'DisplayName', '仿真 P_D');
hold on;
semilogy(SNR_dB, P_D_target_arr, 'b*-', 'DisplayName', '理论 P_D');
xlabel('SNR (dB)');
ylabel('检测概率 P_D');
title(['匹配滤波器检测性能，P_{FA,target} = ', num2str(P_FA_target)]);
legend;
grid on;

% 虚警概率验证图
f2=figure(2);
semilogy(SNR_dB, P_FA_sim, 'k^-', 'DisplayName', '仿真 P_{FA}');
hold on;
yline(P_FA_target, 'r--', 'DisplayName', '目标 P_{FA}');
xlabel('SNR (dB)');
ylabel('虚警概率 P_{FA}');
legend;
title('虚警概率仿真统计');
grid on;
