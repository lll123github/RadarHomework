tau=2;%线性调频信号时间长度
t_all=4;%信号总时间
k=100;%斜率
B=k*tau;%带宽
f_c=300;%中心频率
fs=20*B;%采样率
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
x=A.*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

f1=figure(1);
plot(t,real(x)); % 绘制实部，范围从 -1 到 1
title('x(t)的实部');
xlabel('t');
ylabel('Re(x(t))');

% 线性调频信号的FFT
X = fft(x);
% 绘制FFT幅度谱
f2=figure(2);
plot(f,abs(fftshift(X)),'r'); % 绘制FFT幅度谱
title('x(t)的FFT幅度谱');
xlabel('f');
ylabel('|X(f)|');
% 绘制FFT相位谱
% f3=figure(3);
% plot(f,angle(fftshift(X)),'r'); % 绘制FFT相位谱
% title('x(t)的FFT相位谱');
% xlabel('f');
% ylabel('arg(X(f))');
% 计算匹配滤波器时域信号
match_filter_x=fliplr(conj(x));
y=conv(x,match_filter_x);
% 计算匹配滤波器输出的FFT
Y=fft(y);
f4=figure(4);
f5=figure(5);
f6=figure(6);
f7=figure(7);

for max_phase_error=[0,pi/2,pi]
    % 一次相位误差（线性项）
    a1=max_phase_error/tau;
    x_linear=A.* exp(1i * (pi * k * t.^2+2*pi*f_c*t + a1 * t));

    % 二次相位误差（平方项）
    a2=max_phase_error/tau;
    x_quadratic = A.* exp(1i * (pi * k * t.^2+2*pi*f_c*t + a2 * t.^2));

    % 三次相位误差（立方项）
    a3=max_phase_error/tau;
    x_cubic = A.* exp(1i * (pi * k * t.^2 +2*pi*f_c*t+ a3 * t.^3));

    %随机相位
    rand_phase =max_phase_error * randn(1, N); % 生成随机相位数组
    x_random = A.*exp(1i * (pi * k * t.^2+2*pi*f_c*t + rand_phase));

    %一次相位误差（线性项）
    y_linear=conv(x_linear,match_filter_x);
    y_linear_db=20*log10(abs(y_linear)/max(abs(y_linear)));
    % 二次相位误差（平方项）
    y_quadratic=conv(x_quadratic,match_filter_x);
    y_quadratic_db=20*log10(abs(y_quadratic)/max(abs(y_quadratic)));
    % 三次相位误差
    y_cubic=conv(x_cubic,match_filter_x);
    y_cubic_db=20*log10(abs(y_cubic)/max(abs(y_cubic)));
    % 随机相位
    y_random=conv(x_random,match_filter_x);
    y_random_db=20*log10(abs(y_random)/max(abs(y_random)));

    center_index = floor(length(y)/2);
    range = 500;
    t_index_start = center_index - range;
    t_index_end = center_index + range;

    % 绘制匹配滤波器输出的时域（非DB）
    figure(f4);
    plot(t_new(t_index_start:t_index_end),abs(y_linear(t_index_start:t_index_end)));
    hold on;
    

    figure(f5);
    plot(t_new(t_index_start:t_index_end),abs(y_quadratic(t_index_start:t_index_end)));
    hold on;
    

    figure(f6);
    plot(t_new(t_index_start:t_index_end),abs(y_cubic(t_index_start:t_index_end)));
    hold on;
    

    figure(f7);
    plot(t_new(t_index_start:t_index_end),abs(y_random(t_index_start:t_index_end)));
    hold on;
end

figure(f4);
title('一次相位误差');
xlabel('t');
ylabel('y(t)');
legend('最大相位误差为0','最大相位误差为\pi/2','最大相位误差为\pi');
figure(f5);
title('二次相位误差');
xlabel('t');
ylabel('y(t)');
legend('最大相位误差为0','最大相位误差为\pi/2','最大相位误差为\pi');
figure(f6);
title('三次相位误差');
xlabel('t');
ylabel('y(t)');
legend('最大相位误差为0','最大相位误差为\pi/2','最大相位误差为\pi');
figure(f7);
title('随机相位');
xlabel('t');
ylabel('y(t)');
legend('方差为0','方差为\pi/2','方差为\pi');