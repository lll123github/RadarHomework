tau=2;%线性调频信号时间长度
t_all=4;%信号总时间
k=50;%斜率
B=k*tau;%带宽
% f_c=300;%中心频率
fs=50*B;%采样率
ts=1/fs;%采样间隔
N=floor(t_all*fs)+1;%采样点数(包括端点)

% 定义采样时间点
t= [-(N-1)/2:(N-1)/2]*ts;
% 从 -1 到 1，采样 2001 个点
% 卷积之后的时间点
t_new=[-(N-1):N-1] * ts ;
%频域的绘制范围。频率分辨率应该是fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % 频率范围从 -N/2 到 N/2

%矩形波
A=t<tau/2&t>-tau/2;
% 线性调频信号（无相位误差）
x=A.*exp(1i*(pi*k*t.^2));

I=A.*cos(pi*k*t.^2);
Q=A.*sin(pi*k*t.^2);
S=I+1i.*Q;
match_filter_S=fliplr(conj(S));
y=conv(S,match_filter_S);

center_index = floor(length(y)/2);
range = 100;
t_index_start = center_index - range;
t_index_end = center_index + range;

f1=figure(1);
%计算一般情况下的匹配滤波结果
plot(t_new(t_index_start:t_index_end),abs(y(t_index_start:t_index_end))); % 绘制实部，范围从 -1 到 1

f2=figure(2);
f3=figure(3);
f4=figure(4);
f5=figure(5);
f6=figure(6);
f7=figure(7);
f8=figure(8);
f9=figure(9);
f10=figure(10);

% 幅度不平衡
for epsilon=[0,5,10]
    I_1=A.*(1+epsilon).*cos(pi*k*t.^2);
    Q_1=A.*sin(pi*k*t.^2);
    S_1=I_1+1i.*Q_1;
    y_1=conv(S_1,match_filter_S);
    figure(f2);
    plot(t_new(t_index_start:t_index_end),abs(y_1(t_index_start:t_index_end))); % 绘制实部，范围从 -1 到 1
    hold on;
    figure(f5);
    y_1_db=20*log10(abs(y_1)/max(abs(y_1))); % 转换为dB
    plot(t_new(t_index_start:t_index_end),y_1_db(t_index_start:t_index_end)); % 绘制实部，范围从 -1 到 1
    hold on;
    figure(f8);
    %绘制相位情况
    y_1_phase=angle(y_1);
    plot(t_new(t_index_start:t_index_end),y_1_phase(t_index_start:t_index_end)); % 绘制相位
    hold on;

end
figure(f2);
xlabel('t');
ylabel('|y(t)|');
title('幅度不平衡（电平）');
legend('幅度差为0','幅度差为5','幅度差为10');
figure(f5);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('幅度不平衡（dB）');
legend('幅度差为0','幅度差为5','幅度差为10');
figure(f8);
xlabel('t');
ylabel('phase(y(t))');
title('幅度不平衡（相位）');
legend('幅度差为0','幅度差为5','幅度差为10');


%相位不正交
for delta_theta=[0,pi/4,pi/2]
    I_2=A.*cos(pi*k*t.^2-delta_theta);
    Q_2=A.*sin(pi*k*t.^2);
    S_2=I_2+1i.*Q_2;
    y_2=conv(S_2,match_filter_S);
    figure(f3);
    plot(t_new(t_index_start:t_index_end),abs(y_2(t_index_start:t_index_end))); 
    % 绘制实部，范围从 -1 到 1
    hold on;
    figure(f6);
    y_2_db=20*log10(abs(y_2)/max(abs(y_2))); % 转换为dB
    plot(t_new(t_index_start:t_index_end),y_2_db(t_index_start:t_index_end)); % 绘制实部，范围从 -1 到 1
    hold on;
    figure(f9);
    %绘制相位情况
    y_2_phase=angle(y_2);
    plot(t_new(t_index_start:t_index_end),y_2_phase(t_index_start:t_index_end)); % 绘制相位
    hold on;
end
figure(f3);
xlabel('t');
ylabel('|y(t)|');
title('相位不正交');
legend('相位误差为0','相位误差为\pi/2','相位误差为\pi');
figure(f6);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('相位不正交（dB）');
legend('相位误差为0','相位误差为\pi/4','相位误差为\pi/2');
figure(f9);
xlabel('t');
ylabel('phase(y(t))');
title('相位不正交（相位）');
legend('相位误差为0','相位误差为\pi/4','相位误差为\pi/2');


% 直流偏置

for delta_k=[0,5,10]
    I_3=A.*(1).*cos(pi*k*t.^2)+delta_k;
    Q_3=A.*sin(pi*k*t.^2);
    S_3=I_3+1i.*Q_3;
    y_3=conv(S_3,match_filter_S);
    figure(f4);
    plot(t_new(t_index_start:t_index_end),abs(y_3(t_index_start:t_index_end))); 
    hold on;
    figure(f7);
    y_3_db=20*log10(abs(y_3)/max(abs(y_3))); % 转换为dB
    plot(t_new(t_index_start:t_index_end),y_3_db(t_index_start:t_index_end)); % 绘制实部，范围从 -1 到 1
    hold on;
    figure(f10);
    %绘制相位情况
    y_3_phase=angle(y_3);
    plot(t_new(t_index_start:t_index_end),y_3_phase(t_index_start:t_index_end)); % 绘制相位
    hold on;
end
figure(f4);
xlabel('t');
ylabel('|y(t)|');
title('直流偏置');
legend('同相分量直流偏置为0','同相分量直流偏置为5','同相分量直流偏置为10');

figure(f7);
xlabel('t');
ylabel('20*lg(abs(y(t)))');
title('直流偏置（dB）');
legend('同相分量直流偏置为0','同相分量直流偏置为5','同相分量直流偏置为10');

figure(f10);
xlabel('t');
ylabel('phase(y(t))');
title('直流偏置（相位）');
legend('同相分量直流偏置为0','同相分量直流偏置为5','同相分量直流偏置为10');







