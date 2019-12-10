clc;
clear all;
%% Evolution equation
% i(dA_s/dz)-beta2/2*(d^2A_s/dz^2)+2*gama*abs(A_pump)^2*A_s=0;
%% Time windows
global Tmax Tmin delta_t;
Tmin=-5e-10; Tmax=5e-10; delta_t=1e-13; t=Tmin:delta_t:Tmax;
Ns=floor((Tmax-Tmin)/delta_t+1); df=1/(Tmax-Tmin);  k=-(Ns-1)/2:(Ns-1)/2;  f=fftshift(-k*df);
%% Parameters setting
gama=20e-3;                                 % Nonlinear coefficient (unit: /W*m)
beta2=-10e-27;                              % Group velocity dispersion (unit: s^2/m)
P_pump=1;                                   % Pump peak power (unit: W)
t_pump=16e-12;                              % Half duration of pmup pulse(unit: ps)
%% Characteristic parameters
alfa=sqrt(-4*beta2*gama*P_pump/t_pump^2);   
sigma_t=sqrt(-beta2/alfa);                  % Scale parameter  (unit: s)
T_to_F=sigma_t^2;                           % Time-to-frequency convesion factor (unit: s^2)
Z_T=2*pi/alfa;                              % Self-imaging period (unit: m)
%% Eigen function
% c_l=(2^l*factorial(n)*sqrt(pi))^(-1/2);   % 归一化系数
% pesai_l=c_l*H_l(t/sigma_t)*exp(-t.^2/(2*sigma_t^2));
%% Input field
l=linspace(0,5,6);
hermite0=1;
hermite1=2*t/sigma_t;
hermite2=4*(t/sigma_t).^2-2;
hermite3=8*(t/sigma_t).^3-12*(t/sigma_t);
hermite4=16*(t/sigma_t).^4-48*(t/sigma_t).^2+12;
hermite5=32*(t/sigma_t).^5-160*(t/sigma_t).^3+120*(t/sigma_t);
hermite={'hermite0' 'hermite1' 'hermite2' 'hermite3' 'hermite4' 'hermite5'};

A_input=0;
% Amplitude coefficinents in Fig.2 
a_l=[0.02+0.71i,-0.12-0.25i,0.34+0.56i,-0.12,-0.52+0.81i,0.1+0.58i];
for ii=1:length(l)
    ii
    l_i=l(ii);
    c_l(ii)=(2^l_i*factorial(l_i)*sqrt(pi))^(-1/2);
    
% Generating arbitrary complex amplitude coefficients     
%     r_l(ii)=rand([1,1]);                           
%     phi(ii)=rand([1,1]);
%     a_l(ii)=r_l(i)*exp(1i*phi(i)*2*pi);            

    A_input=A_input+a_l(ii)*c_l(ii)*eval(hermite{ii}).*exp(-t.^2/(2*sigma_t^2));
end
A0=0.2;                                        % Input peak amplitude
A_input=A0*A_input;
%% Input temporal waveform
figure
I_input=abs(A_input).^2;
I_inputmax=max(max(abs(I_input)));
plot(t/1e-12,I_input/I_inputmax,'r-');
xlabel('Time (ps)');ylabel('Normalized Intensity');
xlim([-20,20]);ylim([0,1.2]);
title('Input temporal waveform');
%% Input spectrum envelope
fA_input=fft(A_input);
IfA_input=abs(fA_input).^2;
IfA_inputmax=max(max(abs(IfA_input)));
figure
plot(f/1e9,abs(IfA_input)/IfA_inputmax);
xlabel('Frequency (GHz)');ylabel('Normalized Intensity');
xlim([-500,500]);ylim([0,1.2]);
title('Input spectrum envelope');
%% Evolution in GRIN time lens
Z=8e3;                                      % Paopagation distance (unit: m)                         
num=4e2;                                    % Simulation step number 
dZ=Z/num;                                   % Step length                          
ufft=fft(A_input);
for mm=2:num+1
    mm
  % Linear propagation dZ/2
  linear_operator=1i*beta2*(2*pi*f).^2/2;    % Dispersion opertor
  Dispersion=exp(linear_operator*dZ/2);     
  A_half=ifft(Dispersion.*ufft);
  % Nonlinear propagation dZ
  I_pump=1-t.^2/(t_pump^2);                  % Pump pulse
  XPM_opertor=1i*2*gama*P_pump*I_pump;       % XPM opertor
%   %% Raman_effect
%   t_R=3e-15;                                 % decay time of the Raman gain (unit:s)
%   n_Raman=diff(abs(A_half).^2)/delta_t;
%   n_Raman=[n_Raman,0];
%   T_Raman=exp(-1i*gama*t_R*n_Raman*dZ);
%   %% Self_steepending
%   w0=2*pi*193.55*1e12;                       % center frequency (unit: s^-1)
%   n_SS=diff(A_half.*abs(A_half).^2)/delta_t;
%   n_SS=[n_SS,0];
%   T_SS=exp(-gama/w0*n_SS./A_half*dZ);
%   u1=A_half.*exp(XPM_opertor*dZ).*rectpuls((t)/(2*t_pump)).*T_Raman.*T_SS;
  u1=A_half.*exp(XPM_opertor*dZ).*rectpuls((t)/(2*t_pump)); 
  % linear propagation dZ/2
  ufft=Dispersion.*fft(u1);
  A_out(mm,:)=ifft(ufft);
  I_out(mm,:)=abs(A_out(mm,:)).^2;
  fA_out(mm,:)=ufft;
  IfA_out(mm,:)=abs(fA_out(mm,:)).^2;
end
A_out(1,:)=A_input;
I_out(1,:)=abs(A_input).^2;
fA_out(1,:)=fA_input;
IfA_out(1,:)=abs(fA_input).^2;
%% Temporal intensity patterns
z=0:dZ:Z;
I_outmax=max(max(I_out));
figure
grid off;
mesh(t/1e-12,z/1e3,I_out/I_outmax);    %% 时域分布
view(0,90);
xlim([-15,15]);ylim([0,8]);
xlabel('Time(ps)');ylabel('Propagation distance(km)');
title('Temporal evolution');
colorbar;
%% Spectral intensity patterns
IfA_outmax=max(max(IfA_out));
figure
grid off;
mesh(f/1e9,z/1e3,(IfA_out)/IfA_outmax);
view(0,90);
xlim([-400,400]);ylim([0,8]);
xlabel('Frequency(GHz)'); ylabel('Propagation distance(km)');
title('Spectral evolution');
colorbar;
%% Output time waveform and scaled input spectrum envelope
figure
I0=I_out(45,:);
plot(t/1e-12,abs(I0)/max(max(abs(I0))),'b-');
xlabel('Time (ps)');ylabel('Normalized Intensity');
hold on
plot(2*pi*f*T_to_F/1e-12,abs(IfA_input)/IfA_inputmax,'r-');
xlim([-15,15]);ylim([0,1.4]);
%% 不同位置处频域波形
figure
If0=IfA_out(45,:);
plot(f/1e9,abs(If0)/max(max(abs(If0))),'b-');
xlabel('Frequency (GHz)');ylabel('Normalized Intensity');
hold on
plot(t/T_to_F/(2*pi)/1e9,I_input/I_inputmax,'r-');
xlim([-400,400]);ylim([0,1.4]);