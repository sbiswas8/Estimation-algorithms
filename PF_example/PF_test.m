clc
clear all
close all
profile off
%rng shuffle
w1 = 1e-6;
w2 = 1e-6;
w3 = 1e-6;
Filter_parameter.PF.P_ini = [0.1 0 0; 0 0.01 0; 0 0 1e-4];
Filter_parameter.PF.Q = [w1 0 0; 0 w2 0; 0 0 w3];
Filter_parameter.PF.R = 100^2;
Filter_parameter.Ns = 200;

w1 = 1e-30;
w2 = 1e-30;
w3 = 1e-30;
Filter_parameter.ESP.P_ini = [0.1 0 0; 0 0.01 0; 0 0 1e-4];
Filter_parameter.ESP.Q = [w1 0 0; 0 w2 0; 0 0 w3];
Filter_parameter.ESP.R = 100^2;

%save Filter_parameter Filter_parameter
%run reentry_ref
[PFdata, pfmsgl] = re_entry_PF(Filter_parameter);
[ESPPFdata, espmsgl] = re_entry_ESP_PF(Filter_parameter);
save PFdata PFdata
save ESPPFdata ESPPFdata
t = 0:60;
plot(t,abs(PFdata.error(1:length(t),1))*0.3048,'--k',t,abs(ESPPFdata.error(1:length(t),1))*0.3048,'k')
xlabel('Time (s)')
ylabel('Altitude error (m)')
legend('PF','ESP-PF')
