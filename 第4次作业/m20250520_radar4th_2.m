% ��������
tgtmodel = 'Swerling2';
f = 2e9; % S����Ƶ��
lambda = physconst('LightSpeed')/f; % ��������

% ���߲���
theta_B_deg = 1.1; % ��λ������ȣ��ȣ�
phi_B_deg = 5.7;   % ����������ȣ��ȣ�
theta_B = deg2rad(theta_B_deg);
phi_B = deg2rad(phi_B_deg);
rau = 0.7; % ����Ч��
G_D = (4*pi)/(theta_B*phi_B); % ������
G = G_D * rau; % �������棨����ֵ��
G_dB = 10*log10(G); % �������棨dB��

% �״�ϵͳ����
p_t = 40e3; % ���书�ʣ�W��
tau = 1e-6; % �����ȣ�s��
PRF = 1e3;  % �����ظ�Ƶ�ʣ�Hz��
scan_rate_rpm = 6; % ɨ�����ʣ�ת/���ӣ�
scan_rate = 360*scan_rate_rpm/60; % ɨ�����ʣ���/�룩
T_dwell = theta_B_deg/scan_rate; % פ��ʱ�䣨�룩
n_p = round(PRF * T_dwell); % ����������

% Ŀ�����
RCS = 6; % Ŀ��RCS��m?��
SNR_min_dB = 13; % ���������СSNR��dB��
SNR_min = 10^(SNR_min_dB/10); % ����ֵ

% ϵͳ���
L_s_dB = 3; % ϵͳ��ģ�dB��
L_s = 10^(L_s_dB/10);
L_a_dB = 0; % ������ģ�dB��
L_a = 10^(L_a_dB/10);

% ��������
F_n_dB = 3.5; % ���ջ�����ϵ����dB��
F_n = 10^(F_n_dB/10);
T0 = 290; % �����¶ȣ�K��
k = 1.38e-23; % ������������

% �״﷽�̼������̽�����
numerator = p_t * G^2 * lambda^2 * RCS * n_p * tau;
denominator = (4*pi)^3 * k * T0 * F_n * L_s * L_a * SNR_min;
R_max = (numerator/denominator)^(1/4);

% ��ʾ���
disp(['���̽�����: ' num2str(R_max/1e3) ' km']);

% ʹ��Phased Toolbox��֤
waveform = phased.RectangularWaveform('PulseWidth',tau,'PRF',PRF);
transmitter = phased.Transmitter('PeakPower',p_t,'Gain',G);
target = phased.RadarTarget('MeanRCS',RCS,'OperatingFrequency',f,'Model',tgtmodel);

% ��������źŹ���
[~, range_grid] = rangeaziwaveform(waveform);
tx_sig = transmitter(waveform());
tx_sig = radiator(tx_sig, [0;0]);    % ��������Ŀ�귽��
rx_sig = target(tx_sig);             % Ŀ�귴��
rx_sig = receiver(rx_sig);           % ���ջ�����

% ����ʵ��SNR
noise_bw = 1/tau;
SNR_actual = pow2db(bandpower(rx_sig)/(10*log10(k*Ts*noise_bw) + Fn + L_total));
disp(['ʵ��SNR: ' num2str(SNR_actual) ' dB']);


% ����ͼ���ӻ�����ѡ��
antenna = phased.CosineAntennaElement('FrequencyRange',[f-0.1e9 f+0.1e9],...
    'CosinePower',[1.5 1.5]);
pattern(antenna,f,'Type','powerdb');