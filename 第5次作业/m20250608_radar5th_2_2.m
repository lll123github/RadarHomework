PRF=1e4; % �����ظ�Ƶ�ʣ�Hz��

%�������Ե�Ƶ�ź�
c=3e8; % ���٣���/�룩
tau=2e-8;%���Ե�Ƶ�ź�ʱ�䳤��
t_all=4e-8;%�ź���ʱ��
k=1e11;%б��
B=k*tau;%����
lambda=3e-2;
f_c=c/lambda;%����Ƶ��
% f_c=3e4;%����Ƶ��
fs=20*max(f_c,B);%������
ts=1/fs;%�������
N=floor(t_all*fs)+1;%��������(�����˵�)

% �������ʱ���
t_sample= [-(N-1)/2:(N-1)/2]*ts;
% �� -1 �� 1������ 2001 ����
% ���֮���ʱ���
t_new = [-(N-1):N-1] * ts ;
%Ƶ��Ļ��Ʒ�Χ��Ƶ�ʷֱ���Ӧ����fs/N
f=[-(N-1)/2:(N-1)/2]*fs/N; % Ƶ�ʷ�Χ�� -N/2 �� N/2

%Ŀ��
v=30; % Ŀ���ٶȣ���/�룩
f_d=2*v/lambda; % ������Ƶ��

syms t; % ������ű��� t
% ���δ�
A = @(t) (t > -tau/2) & (t < tau/2);

% �����źţ��ο���
x_tx = @(t) A(t).*exp(1i*(pi*k*t.^2 + 2*pi*f_c*t));

% �����ز���������Ƶƫ + �ӳ�
x_rx = @(t, n) A(t).*exp(1i*(pi*k*t.^2 + 2*pi*(f_c + f_d)*t))+0.1 * randn(size(t));

% ƥ���˲���ģ�壨�ο��źŹ���ʱ�䷴ת��
ref = conj(fliplr(x_tx(t_sample)));
f1=figure(1);
x_sample = x_tx(t_sample); % �����ź�
% plot(t_sample, abs(x_sample), 'DisplayName', 'ʵ��');
% hold on;
% title('���Ե�Ƶ�ź�');
% xlabel('ʱ�� (��)');
% ylabel('����');
% legend('show');


t_accumulate = 0.1;             % ����ʱ�䣨�룩
num_pulses = floor(PRF * t_accumulate);  % �������
compressed_pulses = zeros(num_pulses, 2*N - 1); % ����ѹ���������
for n = 0:num_pulses-1
    % ���ɵ�n���ز�
    echo = x_rx(t_sample);
    % ƥ���˲���������ѹ����
    compressed = conv(echo, ref);
    compressed_pulses(n+1, :) = compressed;
end
% ������ѹ����Ĳ�����ѡ��һ��ʱ��㣨ͨ��ѡ�������Ӧ�ĵ㣩
[~, peak_idx] = max(mean(abs(compressed_pulses), 1));
doppler_data = compressed_pulses(:, peak_idx);  % �ؿ�ʱ������ʱ��FFT

% ������FFT
N_doppler = 2^nextpow2(num_pulses); % 2���ݳ���
doppler_spectrum = fftshift(fft(doppler_data, N_doppler));

% ������Ƶ����
f_d_axis = linspace(-PRF/2, PRF/2, N_doppler);
velocity_axis = f_d_axis * lambda / 2; % �ɶ�����Ƶ��תΪ�ٶ�

% ƥ���˲���άͼ��
figure;
imagesc((1:num_pulses)/PRF, t_new*1e6, 20*log10(abs(compressed_pulses.')));
xlabel('ʱ�� (��)');
ylabel('�ӳ�ʱ�� (��s)');
title('ƥ���˲����������ѹ�������');
colorbar;

% ��������ͼ
figure;
plot(velocity_axis, 20*log10(abs(doppler_spectrum)/max(abs(doppler_spectrum))));
xlabel('�ٶ� (m/s)');
ylabel('���� (dB)');
title('�������ף��ٶȹ��ƣ�');
grid on;


% % ����ѹ��v
% y=conv(x(t_sample),fliplr(conj(x(t_sample))));
% y_max= max(abs(y)); % ��һ��Ŀ������ֵ
% % 3dB ��
% y_3db = 20*log10(0.707); % 3dB��
% % ���ƽ��
% f2=figure(2);
% title('ƥ���˲������');
% xlabel('ʱ�� (��)');
% ylabel('����');
% plot(t_sample,y_3db*ones(size(t_sample)), 'k--', 'DisplayName', '3dB��');
% hold on;
% % ����ƥ���˲������
% plot(t_new, 20*log10(abs(y)./y_max), 'DisplayName', 'ƥ���˲������');
% legend('show');
% hold off;
% % ����Ƶ��
% f3=figure(3);
% title('ƥ���˲��������Ƶ��');
% xlabel('Ƶ�� (Hz)');
% ylabel('������');
% % ����Ƶ��
% Y=fftshift(fft(y,N)); % Ƶ��
% plot(f,abs(Y));
% legend('show');
% hold off;