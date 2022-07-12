function c = CalcPWSCoefficient(data,window_T, dt)
%   Summary of this function goes here.
%   c = CalcPWSCoefficient(data,window_T, dt)
%   Detailed explanation goes here.
%   The function is for calculating the phase-weighted-stack coefficient of N traces.
%   The function is generally adopted in processing flow in Ambient Noise Interferometry
%   for suppressing the incoherent noise.
%
%   IN      
%   data: raw record.
%   window_T: the length of time window (s).
%   dt: the interval (s) in time domain of recorded data.
%
%  OUT   
%   c: output result, the phase-weighted-stack coefficient, it's used to ajust
%      the instantaneous amplitude of stacked trace.
%
%  References:
%  Schimmel, M., & Paulssen, H. (1997). Noise reduction and detection of weak, 
%  coherent signals through phase weighted stacks. Geophysical Journal International, 
%  130, 497¨C505. https://doi.org/10.1111/j.1365-246X.1997.tb05664.x
%
%  Author(s): Yan Yingwei
%  Copyright: 2020-2025 
%  Revision: 1.0  Date: 7/27/2020
%
%  Academy of Opto-Electronics, China Electronic Technology Group Corporation (AOE CETC)

[M,N] = size(data);
P = floor(0.5*window_T/dt);
c = zeros(M,1);

data_hilbert = hilbert(data);
data_hilbert_angle = angle(data_hilbert);

for i=1:M
    c(i) = abs(sum(exp(1i*data_hilbert_angle(i,:))));
end
c = c/N;
c_temp = c;

for i=1:P
    P_temp = (i-1);
    c(i) = sum(c_temp(i-P_temp:i+P_temp))/(2*P_temp+1);
end

for i=P+1:M-P
    c(i) = sum(c_temp(i-P:i+P))/(2*P+1);
end

for i=M-P+1:M
    P_temp = M-i;
    c(i) = sum(c_temp(i-P_temp:i+P_temp))/(2*P_temp+1);
end
end

