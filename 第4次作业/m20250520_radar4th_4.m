% ��������
tgtmodel = 'Swerling2';
f = 2e9; % S����Ƶ��
c= physconst('LightSpeed'); % ����
lambda = c/f; % ��������

% ���߲���
theta_B_deg = 1.1; % ��λ������ȣ��ȣ�
phi_B_deg = 5.7;   % ����������ȣ��ȣ�
theta_B = deg2rad(theta_B_deg);
phi_B = deg2rad(phi_B_deg);
rau = 0.7; % ����Ч��
G_D_max = (4*pi)/(theta_B*phi_B); % ������
G_max = G_D_max * rau; % �������棨����ֵ��
G_dB_max = 10*log10(G_max); % �������棨dB��
arrayLength = [7.3 1.52]; % ���г��ȣ��ף�

%���߽�ģ
%������Ԫ���
d=sqrt(arrayLength(1)*arrayLength(2)*theta_B*phi_B/(4*pi)); % ��������Ϊ��������
%����d��lamda�ı�ֵ
d_ratio = d/lambda; % ��Ԫ����벨���ı�ֵ��
% ��ӡ���
disp(['��Ԫ���: ' num2str(d) ' ��']);
disp(['��Ԫ����벨���ı�ֵ: ' num2str(d_ratio)]);

% ������Ԫ�������������ھ��Ⱦ�������ģ�ͣ�
N_az = round(arrayLength(1)/d); % ��λ����Ԫ��
N_el = round(arrayLength(2)/d); % ��������Ԫ��

G=arrayGain(theta, phi, N_az, N_el, lambda, d, G_max);
f1=figure(1);
imagesc(rad2deg(theta), rad2deg(phi), 10*log10(G));
xlabel('Azimuth Angle (deg)');
ylabel('Elevation Angle (deg)');
title('2D Gain Pattern (dB)');
colorbar;
axis xy;

% �״�ϵͳ����
p_t = 40e3; % ���书�ʣ�W��
tau = 1e-6; % �����ȣ�s��

scan_rate_rpm = 6; % ɨ�����ʣ�ת/���ӣ�
scan_rate = 360*scan_rate_rpm/60; % ɨ�����ʣ���/�룩
T_dwell = theta_B_deg/scan_rate; % פ��ʱ�䣨�룩


% Ŀ�����
RCS = 6; % Ŀ��RCS
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

denominator = (4*pi)^3 * k * T0 * F_n * L_s * L_a * SNR_min;


PRF_list=[1e4,2e4, 4e4, 1e5]; % �����ظ�Ƶ���б�
R_max_list = zeros(size(PRF_list)); % ��ʼ�����̽������б�
fig_index=2;
for i = 1:length(PRF_list) 
    PRF = PRF_list(i);
    n_p = round(PRF * T_dwell); % ����������
    % �״﷽�̼������̽�����
    numerator = p_t * G_max^2 *( lambda^2 * RCS * n_p * tau);
    R_max = (numerator/denominator)^(1/4);
    % ��ʾ���
    disp(['�����ظ�Ƶ��: ' num2str(PRF/1e3) ' kHz, ���̽�����: ' num2str(R_max/1e3) ' km']);
    R_max_list(i) = R_max; % �洢���
    %��ͼ
    figure(fig_index);
    plot(PRF_list/1e3, R_max_list/1e3, '-o');
    xlabel('�����ظ�Ƶ�� (kHz)');
    ylabel('���̽����� (km)');
    title(['PRF = ', num2str(PRF/1e3), ' kHz']);
end
% ��ͼ(����)
f2=figure(2);
plot(PRF_list/1e3, R_max_list/1e3, '-o');
set(gca, 'XScale', 'log'); % ����x��Ϊ�����̶�
set(gca, 'YScale', 'linear'); % ����y��Ϊ���Կ̶�
xlabel('�����ظ�Ƶ�� (kHz)');
ylabel('���̽����� (km)');




% PRF = 1e5;  % �����ظ�Ƶ�ʣ�Hz����Ҫ���ǲ�ͬ�����Ľ��
% n_p = round(PRF * T_dwell); % ����������
% % �״﷽�̼������̽�����
% numerator = p_t * G^2 * lambda^2 * RCS * n_p * tau;
% R_max = (numerator/denominator)^(1/4);
% % ��ʾ���
% disp(['���̽�����: ' num2str(R_max/1e3) ' km']);




% % ����ͼ���ӻ�����ѡ��
% antenna = phased.CosineAntennaElement('FrequencyRange',[f-0.1e9 f+0.1e9],...
%     'CosinePower',[1.5 1.5]);
% pattern(antenna,f,'Type','powerdb');

% % ����theta��phi����
% [theta, phi] = meshgrid(linspace(-pi/6, pi/6, 100), linspace(-pi/2, pi/2, 100));

% % ��������ͼģ�ͣ������˹��cos�ͣ�
% gain = cos(theta/(theta_B/2)).^2 .* cos(phi/(phi_B/2)).^2;
% gain(gain < 0) = 0;
% gain_dB = 10*log10(gain);
% gain_dB(gain == 0) = -60;

% % ��ͼ
% figure;
% surf(rad2deg(theta), rad2deg(phi), gain_dB, 'EdgeColor', 'none');
% xlabel('Azimuth (deg)');
% ylabel('Elevation (deg)');
% zlabel('Gain (dB)');
% title('2D Radiation Pattern');
% colorbar;


function G = arrayGain(theta, phi, N_az, N_el, lambda, d, G_max)
        % theta, phi: �Ƕ� (rad)
        % N_az, N_el: ��λ�븩��������Ԫ��
        % lambda: ����
        % d: ��Ԫ���
        % G_max: ������棨����ֵ��

        k = 2*pi/lambda;

        % ��λ��
        psi_az = k * d * sin(theta) .* cos(phi);
        psi_el = k * d * sin(phi);

        % ��ֹ����
        psi_az(psi_az==0) = 1e-12;
        psi_el(psi_el==0) = 1e-12;

        % ��һ����������
        AF_az = sin(N_az * psi_az / 2) ./ (N_az * sin(psi_az / 2));
        AF_el = sin(N_el * psi_el / 2) ./ (N_el * sin(psi_el / 2));

        AF = AF_az .* AF_el;

        % ���溯��
        G = G_max * abs(AF).^2;
end