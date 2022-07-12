function [E,freq,v] = CLPhaseShiftofSW(rec,dt,r,RangeFreq,ParaVelocity,mId,tensorName,dataMode)
%   Summary of this function goes here.
%   [E,freq,v] = CLPhaseShiftofSW(rec,dt,r,RangeFreq,ParaVelocity,mId,tensorName,dataMode)
%   Detailed explanation goes here.
%   The function is for generating the dispersive energy of surface waves
%   based on the cylindrical wave diffusion.
%
%   IN      
%           rec: the cross-correlation record of surface waves, the size of 
%                rec is [nt*nx].
%            dt: the sampling interval in time domain (s).
%             r: it's a vector, the distance (m) of the receiver/station point 
%                to the reference point.
%     RangeFreq: the range of frequency (Hz) in the dispersive energy, such
%                as [5 100].
%  ParaVelocity: the parameter of the phase velocity (m/s), first is the 
%                minimal, second is the maximal, third is the interval, such 
%                as [50 800 1].
%           mId: the identification for the phase velocity measurement,
%                'basic-function' means the measurement based on the sine 
%                and cosine basis function,'magnitude-based' means
%                cylindrical-wave phase.
%    tensorName: 'ZZ','RR'.
%      dataMode: '0' reprensents the common multi-record of surface waves,in 
%                this case, the broadband source is satisfied; '1' represents
%                the data excited by the rammer.
%                
%
%  OUT   
%            E:  the normalized dispersive energy.
%          freq: the frequency (Hz) vecotor of the normalized dispersive energy.
%             v: the phase velocity (m/s) vector of the normalized dispersive energy.
%  
%  References: 
%  Zywicki, D., & Rix, G. 2005. Mitigation of the Near-Field Effects for 
%  Seismic Surface Wave Velocity Estimation with Cylindrical Beamformers, 
%  Journal of Geotechnical and Geoenviromental Engineering, 131(8), 970-977,
%  https://doi.org/10.1061/(ASCE)1090-0241(2005)131:8(970).
%
% Haney, M., & Nakahara, H. 2014. Surface-Wave Green¡¯s Tensors in the 
% Near Field, Bulletin of the Seismological Society of America, 104(3), 
% 1578-1586, https://doi.org/10.1785/0120130113.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2022-2025 
%  Revision: 1.0  Date: 4/8/2022
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

if nargin==7
    dataMode = 0;
elseif nargin==6
    dataMode = 0;
    tensorName = 'ZZ';
elseif nargin==5
    dataMode = 0;
    tensorName = 'ZZ';
    mId = 'basic-function';
end

fmin = RangeFreq(1);
fmax = RangeFreq(2);
vmin = ParaVelocity(1);
vmax = ParaVelocity(2);
dv = ParaVelocity(3);

[nt,nx] = size(rec);       % the shape of CCFs
rec_fft = zeros(nt,nx);    % spectrums 

v = vmax:-dv:vmin;         % scanned phase velocity vector
lv = length(v);            % the length of variable 'v'

Fs = 1/dt;                 % sampling rate, Hz
f = Fs*(0:(nt/2))/nt;      % frequency vector

% the index of the minimum frequency of the dispersion image
ind = find(f>fmin);       
fmin_ind  = ind(1)-1;

% the index of the maximum frequency of the dispersion image
ind = find(f>fmax);
fmax_ind = ind(1);

freq = f(fmin_ind:fmax_ind); % frequency vector of the dispersion image
lf = length(freq);

for j=1:nx
    rec_fft(:,j) = fft(rec(:,j),nt,1); % get the spectrums of the CCFs
end


rec_fft_amp = abs(rec_fft);
rec_fft_n = rec_fft./rec_fft_amp; % normalized the spectrums of the CCFs

rec_fft_n_r = real(rec_fft_n); % real part of the normalized spectrums
rec_fft_n_i = imag(rec_fft_n); % imaginary part of the normalized spectrums

E = zeros(lv,lf); % the matrix of the dispersive energy

% calcualting the dispersive energy of the CCFs
if strcmp(tensorName,'ZZ')
    if strcmp(mId,'basic-function')
        for j=1:lf
            f_ind = fmin_ind+j-1;
            for i=1:lv
                ko = 2*pi*f(f_ind)/v(i);
                for k=1:nx
                    x = r(k);
                    H0 = besselh(0,2,ko*x);
                    realPart = real(H0);
                    imagPart = imag(H0);
                    phi = atan2(imagPart,realPart);
                    %                 phi = -ko*x+pi/4;
                    E(i,j) = E(i,j)+cos(phi)*rec_fft_n_r(f_ind,k)+sin(phi)*rec_fft_n_i(f_ind,k);
                end
            end
        end
    E(E<0.0) = 0.0;
    elseif strcmp(mId,'magnitude-based')
        for j=1:lf
            f_ind = fmin_ind+j-1;
            for i=1:lv
                ko = 2*pi*f(f_ind)/v(i);
                for k=1:nx
                    x = r(k);
                    H0 = besselh(0,2,ko*x);
                    realPart = real(H0);
                    imagPart = imag(H0);
                    phi = atan2(imagPart,realPart);
                    E(i,j) = E(i,j)+exp(-1i*phi)*rec_fft_n(f_ind,k);
                end
            end
        end
        E = abs(E);
    end
end

if strcmp(tensorName,'RR')
    if strcmp(mId,'basic-function')
        for j=1:lf
            f_ind = fmin_ind+j-1;
            for i=1:lv
                ko = 2*pi*f(f_ind)/v(i);
                for k=1:nx
                    x = r(k);
                    H02 = besselh(0,2,ko*x);
                    H22 = besselh(2,2,ko*x);
                    HH = H02-H22;
                    realPart = real(HH);
                    imagPart = imag(HH);
                    phi = atan2(imagPart,realPart);
%                     phi = -ko*x;
                    E(i,j) = E(i,j)+cos(phi)*rec_fft_n_r(f_ind,k)+sin(phi)*rec_fft_n_i(f_ind,k);
                end
            end
        end
    E(E<0.0) = 0.0;
    elseif strcmp(mId,'magnitude-based')
        for j=1:lf
            f_ind = fmin_ind+j-1;
            for i=1:lv
                ko = 2*pi*f(f_ind)/v(i);
                for k=1:nx
                    x = r(k);
                    H02 = besselh(0,2,ko*x);
                    H22 = besselh(2,2,ko*x);
                    HH = H02-H22;
                    realPart = real(HH);
                    imagPart = imag(HH);
                    phi = atan2(imagPart,realPart);
%                     phi = -ko*x+5*pi/4;
                    E(i,j) = E(i,j)+exp(-1i*phi)*rec_fft_n(f_ind,k);
%                      E(i,j) = E(i,j)+(exp(-1i*phi1)-exp(-1i*phi2))*rec_fft_n(f_ind,k);
                end
            end
        end
        E = abs(E);
    end
end

if dataMode==0
    for j=1:lf
        E(:,j) = E(:,j)./max(E(:,j));
    end
else
    E = E./max(max(E));
end
E(isnan(E)) = 0.0;
E(isinf(E)) = 0.0;
E = E.*E;
end

