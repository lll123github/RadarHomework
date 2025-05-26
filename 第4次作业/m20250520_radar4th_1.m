%RCS��ģ Swerling2ģ��
tgtmodel = 'Swerling2';
f=2e9; % ����Ƶ�� S����
antenna = phased.IsotropicAntennaElement('BackBaffled',true);
radiator = phased.Radiator('OperatingFrequency',f,'Sensor',antenna);
RCS=6;
target = phased.RadarTarget('MeanRCS',RCS,'OperatingFrequency',f,...
    'Model',tgtmodel);

sigma = target.MeanRCS; % Ŀ���״�ɢ�����

p_t=40*10e3; % ���书��

lambda=1/f;% ����
theta_B=1.1/360*pi;
phi_B=5.7/360*pi; % �������
G_D=4*pi/(theta_B*phi_B);
rau=0.7; % ����Ч��
G=G_D*rau; % ��������

SN_min_db=13; % �����
SN_min=10^(SN_min_db/10); % ��С�����
L_s=10^(L_s_dB/10); % ϵͳ���
% L_a=10^(L_a_dB/10); % �������
L_a=1; % �������

p_r=(p_t*G^2*lambda^2*sigma)/((4*pi)^3*L_s*L_a); % ���չ���

PRF=1e3; % �����ظ�Ƶ��
tau=1*10e-6; % ������
% f_s=1e11; % ����Ƶ��
% B_n=f_s; % ��������
T=6/60; % �����ظ����� ÿ������ת
T_0=290; % �ο��¶�
F_n_dB=3.5; % ����ϵ��(dB)
F_n=10^(F_n_dB/10); % ����ϵ��
k=1.38e-23; % ������������

R_max=(p_r*tau*PRF*T/(F_n*k*T_0*SN_min))^(1/4); % ���̽�����