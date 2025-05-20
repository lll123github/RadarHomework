tau=2;%���Ե�Ƶ�ź�ʱ�䳤��
t_all=4;%�ź���ʱ��
k=100;%б��
B=k*tau;%����
f_c=300;%����Ƶ��
fs=20*B;%������
ts=1/fs;%�������
N=floor(t_all*fs)+1;%��������(�����˵�)

% �������ʱ���
t= [-(N-1)/2:(N-1)/2]*ts;
% �� -1 �� 1������ 2001 ����
% ���֮���ʱ���
t_new=[-(N-1):N-1] * ts ;
%Ƶ��Ļ��Ʒ�Χ��Ƶ�ʷֱ���Ӧ����fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % Ƶ�ʷ�Χ�� -N/2 �� N/2

%���β�
A=t<tau/2&t>-tau/2;
% ���Ե�Ƶ�źţ�����λ��
x=A.*exp(1i*(pi*k*t.^2+2*pi*f_c*t));

f1=figure(1);
plot(t,real(x)); % ����ʵ������Χ�� -1 �� 1
title('x(t)��ʵ��');
xlabel('t');
ylabel('Re(x(t))');

% ���Ե�Ƶ�źŵ�FFT
X = fft(x);
% ����FFT������
f2=figure(2);
plot(f,abs(fftshift(X)),'r'); % ����FFT������
title('x(t)��FFT������');
xlabel('f');
ylabel('|X(f)|');
% ����FFT��λ��
% f3=figure(3);
% plot(f,angle(fftshift(X)),'r'); % ����FFT��λ��
% title('x(t)��FFT��λ��');
% xlabel('f');
% ylabel('arg(X(f))');
% ����ƥ���˲���ʱ���ź�
match_filter_x=fliplr(conj(x));
y=conv(x,match_filter_x);
% ����ƥ���˲��������FFT
Y=fft(y);
f4=figure(4);
f5=figure(5);
f6=figure(6);
f7=figure(7);

for max_phase_error=[0,pi/2,pi]
    % һ����λ�������
    a1=max_phase_error/tau;
    x_linear=A.* exp(1i * (pi * k * t.^2+2*pi*f_c*t + a1 * t));

    % ������λ��ƽ���
    a2=max_phase_error/tau;
    x_quadratic = A.* exp(1i * (pi * k * t.^2+2*pi*f_c*t + a2 * t.^2));

    % ������λ�������
    a3=max_phase_error/tau;
    x_cubic = A.* exp(1i * (pi * k * t.^2 +2*pi*f_c*t+ a3 * t.^3));

    %�����λ
    rand_phase =max_phase_error * randn(1, N); % ���������λ����
    x_random = A.*exp(1i * (pi * k * t.^2+2*pi*f_c*t + rand_phase));

    %һ����λ�������
    y_linear=conv(x_linear,match_filter_x);
    y_linear_db=20*log10(abs(y_linear)/max(abs(y_linear)));
    % ������λ��ƽ���
    y_quadratic=conv(x_quadratic,match_filter_x);
    y_quadratic_db=20*log10(abs(y_quadratic)/max(abs(y_quadratic)));
    % ������λ���
    y_cubic=conv(x_cubic,match_filter_x);
    y_cubic_db=20*log10(abs(y_cubic)/max(abs(y_cubic)));
    % �����λ
    y_random=conv(x_random,match_filter_x);
    y_random_db=20*log10(abs(y_random)/max(abs(y_random)));

    center_index = floor(length(y)/2);
    range = 500;
    t_index_start = center_index - range;
    t_index_end = center_index + range;

    % ����ƥ���˲��������ʱ�򣨷�DB��
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
title('һ����λ���');
xlabel('t');
ylabel('y(t)');
legend('�����λ���Ϊ0','�����λ���Ϊ\pi/2','�����λ���Ϊ\pi');
figure(f5);
title('������λ���');
xlabel('t');
ylabel('y(t)');
legend('�����λ���Ϊ0','�����λ���Ϊ\pi/2','�����λ���Ϊ\pi');
figure(f6);
title('������λ���');
xlabel('t');
ylabel('y(t)');
legend('�����λ���Ϊ0','�����λ���Ϊ\pi/2','�����λ���Ϊ\pi');
figure(f7);
title('�����λ');
xlabel('t');
ylabel('y(t)');
legend('����Ϊ0','����Ϊ\pi/2','����Ϊ\pi');