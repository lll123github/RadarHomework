%产生线性调频信号
c=3e8; % 光速（米/秒）
tau=2e-8;%线性调频信号时间长度
k=2e17;
B=k*tau;%带宽
lambda=3e-3;
f_c=c/lambda;%中心频率
fs=5*max(f_c,B);%采样率

% 计算多普勒频移
R_target = 1e5; % 第一个目标距离（米） %为了规避大数吃小数，把R_target调成0
v=300; % 目标速度（米/秒）
f_d=2*v/lambda;
% fs=4*f_d; % 采样率
ts=1/fs;%采样间隔
PRF=1e4; % 脉冲重复频率（Hz）
t_accumulate=0.1; % 脉冲串累积时间
N_accumulate=floor(t_accumulate*PRF);%叠加信号数(包括端点)

%矩形波
A=@(t) t<tau/2&t>-tau/2;
% 线性调频信号（无相位误差）
x=@(t) A(t).*exp(1i*(pi*k*t.^2+2*pi*f_c*t));
% 快时间回波数
N_fast=10;
N_slow=N_accumulate/N_fast; % 慢时间回波的组数
% 慢时间回波的组数
%计算回波信号：



% 单个回波信号的采样时间点
%最早回来的时刻
t_early=2*R_target/c; % 最早回来的时刻
% 最晚回来的时刻
t_late=2*(R_target+v*t_accumulate)/c; % 最晚回来的时刻
N_s_below=floor(t_early/ts);
N_s_top=floor(t_late/ts);
N_edge=floor((t_late-t_early)/ts/4);
s_t_sample=[N_s_below:N_s_top]*ts; % 从0开始采样

% 不建议使用匿名函数的方法，因为需要对函数进行深拷贝进行额外的定义
% F=cell(1,N_slow);
% F={};
% for i = 0:N_slow-1
%     s=@(t) 0;
%     for j= 0:N_fast-1
%         t_accumulated=(i*N_fast+j)*t_accumulate/N_accumulate; % 已累积时间
%         R_delta=v*t_accumulated; % 目标距离
%         % 计算回波时刻
%         t_ret=2*(R_target+R_delta)/c; % 回波时间
%         % 计算多普勒频移
%         f_d=2*v/lambda;
%         % 多普勒频移
%         x_d=@(t) x(t).*exp(1i*2*pi*f_d*t);
%         %需要注意变量的定义范围，以避免引用的发生
%         delay=t_ret+t_accumulated;
%         s_old=copy(s);
%         % 计算回波信号
%         s =@(t) s_old(t)+x(t-copy(delay));
%     end
%     % 在F中记录回波信号
%     F{end+1} = s;
% end
F=zeros(N_slow,length(s_t_sample));
G=zeros(N_fast,N_slow,length(s_t_sample));
for i=0:N_slow-1
    s=zeros(1,length(s_t_sample));
    for j=0:N_fast-1
        t_accumulated=(i*N_fast+j)*t_accumulate/N_accumulate; % 已累积时间
        R_delta=v*t_accumulated; % 目标距离
        % 计算回波时刻
        t_ret=2*(R_target+R_delta)/c; % 回波时间
        % delay=t_ret+t_accumulated;
        % 计算多普勒频移
        f_d=2*v/lambda;
        % 多普勒频移
        x_d=@(t) x(t-t_ret).*exp(1i*2*pi*f_d*(t-t_ret));
        % 计算回波信号
        delta_s=sum(x_d(s_t_sample));
        s =s+x_d(s_t_sample);% 快时间累积信号
        G(j+1,i+1,:)=x_d(s_t_sample); % 每一个脉冲的回波信号
    end
    F(i+1,:)=s;
end

% 对快时间的脉冲累积做匹配滤波
matched_filter_1= @(t) fliplr(conj(x(t)));

% 定义采样时间点
% 匹配滤波器的采样时间点
t_all=2*tau;
N_x=floor(t_all*fs)+1;%采样点数(包括端点)
N_x_below=-(N_x-1)/2;
N_x_top=(N_x-1)/2;
x_t_sample= [N_x_below:N_x_top]*ts; % 从0开始采样



% 匹配滤波之后的采样时间点
y_t_sample= [N_x_below+N_s_below:N_x_top+N_s_top]*ts;

% 对一组快时间的回波信号进行匹配滤波
s1=F(1,:);
y1 = conv(s1, matched_filter_1(x_t_sample)); % 计算匹配滤波器输出
y1_max = max(abs(y1));
y1_db= 20*log10(abs(y1)/y1_max); % 转换为dB
f1=figure(1);
% 绘制匹配滤波器输出
plot(y_t_sample, abs(y1), 'DisplayName', '匹配滤波器输出');
title('匹配滤波器输出');
xlabel('时间 (秒)');
ylabel('幅度');

% 分别做脉冲压缩
N=zeros(N_slow,length(y_t_sample));
for i=1:N_slow
    s=F(i,:);
    y= conv(s, matched_filter_1(x_t_sample)); % 计算匹配滤波器输出
    N(i,:)=y;
end

% 按列叠加匹配滤波之后的结果
y_total=sum(N,1);
y_total_max = max(abs(y_total));
y_total_db= 20*log10(abs(y_total)/y_total_max); % 转换为dB
f2=figure(2);
% 绘制匹配滤波器输出
plot(y_t_sample, abs(y_total), 'DisplayName', '匹配滤波器输出');
title('匹配滤波器输出');
Q=zeros(N_slow,length(y_t_sample));
% 对N做慢时间的FFT
for j=1:length(y_t_sample)
    Q(:,j)=abs(fftshift(fft(N(:,j))));
end
f3=figure(3);
% 绘制二维幅度谱
y_f_sample=linspace(-fs/2, fs/2, length(y_t_sample)); % 频率采样点
y_t_truth=y_t_sample+2*R_target/c; % 真实的回波时间
imagesc(y_t_truth,y_f_sample,Q);
xlabel('接收延迟时间(秒)');
ylabel('频率 (Hz)');
title('对快时间信号进行叠加后的幅度谱');


% 对G进行处理
G_total=sum(G,2);
% 得到100条结果
% 取第一条进行匹配滤波
G1=G_total(1,1,:);
G1=reshape(G1,1,length(s_t_sample));
% 解调信号
h = @(t) exp(-1i*2*pi*f_c*t);
% 进行解调
G1=G1.*h(s_t_sample);
matched_filter_2= @(t) fliplr(conj(x(t).*h(t))); % 匹配滤波器
% 匹配滤波
Y=conv(G1, matched_filter_2(x_t_sample));
Y_max = max(abs(Y));
Y_db= 20*log10(abs(Y)/Y_max);
f4=figure(4);
plot(y_t_sample, abs(Y), 'DisplayName', '匹配滤波器输出');
title('匹配滤波器输出');
xlabel('时间 (秒)');
ylabel('幅度');


%做FFT
Y_fft=abs(fftshift(fft(Y)));
f5=figure(5);
% 绘制二维幅度谱
y_f_sample=linspace(-fs/2, fs/2, length(y_t_sample));
Y_fft_db=20*log10(Y_fft/max(Y_fft));
plot(y_f_sample,Y_fft_db);
xlabel('频率 (Hz)');
ylabel('幅度');
xlim([0,2*f_d]);


% 得到幅度最高的那个频率
[~,index]=max(Y_fft);
fc_index=floor((f_c/fs+0.5) *length(y_f_sample));
delta_f=abs(index-fc_index);

% 频率分辨率不够用来显示多普勒频移

Y_total=zeros(N_fast,length(y_t_sample));
% for i=1:N_fast
for i=1:1
    % 对每一组慢时间的回波信号进行匹配滤波
    s=G(1,i,:).*exp(-1i*2*pi*f_c*s_t_sample); % 解调
    Y_total(i,:)=Y_total(i,:)+s;
    % 匹配滤波
    Y=conv(s, matched_filter_2(x_t_sample));
end
%解调
Y_fft=abs(fftshift(fft(Y)));
Y_fft_db=20*log10(Y_fft/max(Y_fft));
f6=figure(6);
plot(y_t_sample, Y_fft_db);
xlabel('频率');
ylabel('幅度 (dB)');
title('频谱');