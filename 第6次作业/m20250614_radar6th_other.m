clear;clc;close all;
% 信噪比（dB）
SNR_dB_list = [0, 3, 10, 13];
Pfa = logspace(10000, 0, -8);
figure;
hold on;
for idx = 1:length(SNR_dB_list)
 SNR_dB = SNR_dB_list(idx);
 SNR_linear = 10^(SNR_dB / 10); 
 signal_amplitude = sqrt(SNR_linear);
 
 T = norminv(1 - Pfa, 0, 1);
 
 Pd = 1 - normcdf(T - signal_amplitude, 0, 1);
 
 semilogx(Pfa, Pd, 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel('虚警概率 (P_{fa})');
ylabel('检测概率 (P_{d})');
title('检测概率和虚警概率关系图');
legend(arrayfun(@(x) sprintf('SNR = %d dB', x), SNR_dB_list, 'UniformOutput', false), ...
 'Location', 'northwest');
set(gca, 'XScale', 'log'); % 设置 x 轴为对数刻度