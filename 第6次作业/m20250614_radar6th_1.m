% 高斯白噪声下的目标检测仿真
clc; clear;

% 信噪比范围（单位 dB）
SNR_dB = [0,3,10,13];
SNR_lin = 10.^(SNR_dB / 10);  % 转换为线性 S/N

% 设置门限 T 的范围
T_high = logspace(0, 2, 500);
T_low=linspace(-3, 1, 100);
T=[T_low,T_high]; % 合并两个范围

% 初始化结果存储
P_FA_all = zeros(length(SNR_lin), length(T));
P_D_all = zeros(length(SNR_lin), length(T));

% 计算不同SNR下的 P_FA 与 P_D
for k = 1:length(SNR_lin)
    snr = SNR_lin(k);
    sqrt_snr = sqrt(snr);
    
    % % 虚警概率
    % P_FA = 0.5 * (1 - erf(T ./ sqrt_snr));
    
    % % 检测概率
    % P_D = 0.5 * (1 - erf(T - sqrt_snr));

    % 虚警概率
    P_FA=1-normcdf(T,0,sqrt_snr);

    % 检测概率
    P_D=1-normcdf(T,snr,sqrt_snr);

    %PPT 上面的函数是不正确的，使用的erf和matlab定义的erf是不一样的
    
    % 存储
    P_FA_all(k, :) = P_FA;
    P_D_all(k, :) = P_D;
end

% 绘图：以 10 - log10(P_FA) 为横坐标
figure;
hold on;
colors = lines(length(SNR_dB));

for k = 1:length(SNR_dB)
    x_plot =- log10(P_FA_all(k, :));  % 横坐标转换
    % x_plot=(P_FA_all(k, :));
    plot(x_plot, P_D_all(k, :), 'Color', colors(k,:), ...
        'DisplayName', ['SNR = ', num2str(SNR_dB(k)), ' dB']);
end

xlabel('-log_{10}(P_{FA})');
ylabel('检测概率 P_D');
title('高斯白噪声下的检测性能');
legend show;
grid on;
xlim([0 8]);  
ylim([0 1]);
