%% Multimodal Dispersion Measurement for the Long Line I DAS data
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-
% Purpose of this script:
% This is a script for measuring the multimodal dispersion of the Long Line
% I distributed-acoustic-sensing (DAS) data based on the cylindrical-wave
% phase. And the DAS data sets are collected at the Garner Valley Downhole 
% Array (GVDA) field test site in Southern California in September 2013, 
% and the survey area is adjacent to Highway 74 of California. 
%
% About the information of survey area and data, please visit:
% https://gdr.openei.org/submissions/481
% http://nees.ucsb.edu/facilities/GVDA
%
% It reproduces the dispersion measurement results of the article 2.
%
% Note:
% Before using the script, we strongly recommend that you read our article 
% and the screen recording tutorial. The article and tutorial are placed in
% "main-folder"\Doc. If you publish an article by using our programs, 
% please cite our article. If there are any questions about the code, 
% please contact with the first author Yan Yingwei by email 
% "wallace2012y@outlook.com".
%
% Please ensure the main folder is in the current directory of 
% matlab, and its subfolders are also added the path of matlab.
%
% The code has been debugged on the the matlab R2019a. A higher or lower 
% version of matlab may also make the code run successfully.
%
% It would run for about a few hours.
%
% References: 
% Haney, M., & Nakahara, H. 2014. Surface-Wave Green��s Tensors in the 
% Near Field, Bulletin of the Seismological Society of America, 104(3), 
% 1578-1586, https://doi.org/10.1785/0120130113.
%
% Our article:
% Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
% the multimodal surface wave dispersion curves: Staging strategy and Pattern 
% search with embedded Kuhn-Munkres algorithm, Geophysical Journal
% International,231(01), 47-71, 
% https://doi.org/10.1093/gji/ggac178. 
%
% Yan， Y., Chen, X., Li, J., Guan, J., Xi, C., Liu, H. 2023. Inversion of
% multimodal dispersion curves from distributed acoustic sensing measurements
% for subsurface imaging: A field case of Garner Valley, California, Journal of
% Applied Geophysics, 214, 105070,
% https://doi.org/10.1016/j.jappgeo.2023.105070 
%
% Author(s): Yan Yingwei
% Email:     wallace2012y@outlook.com
% Copyright: 2022-2025 
% Revision:  1.0  Date: 9/6/2023
%
% Department of Earth and Space Sciences, Southern University of Science 
% and Technology (SUSTech).
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-

clear;                  % clear the workspace of matlab
%% load the DAS data of Long Line I
currDir = cd;
filePath = strcat(currDir,'\Data\LongLineIData_GVDA.mat');
load(filePath);

data = LongLineIData_GVDA; % each column represent a channel of DAS data
Fs = 200;            % sampling rate, Hz
dt = 1/Fs;           % sampling interval, s

[M,N] = size(data); % the shape of raw DAS data
T = M*dt;           % the total recording time, s
%% set the parameters for getting all CCFs gathers
secT =5;            % the time-length of CCF��s
alp = 4/5;          % overlap ration of adjacent time-window
secPT = floor(secT/dt);    % the number of samples of each CCF

% the number of segments of the CCF, or stacking 
secNum = floor((T-alp*secT)/(secT*(1-alp)));

% let each sub-window contain 60 DAS channels
virTraceNum = 59;            % each CCFs gather contains 59 CCFs

% total number of CCFs gather
virShotNum = N- virTraceNum; 

CCFsC = zeros(secPT,virTraceNum);   % causal part of CCFs, temp variable
CCFsAC = zeros(secPT,virTraceNum);  % acausal part of CCFs, temp variable

% the causal part of all the CCFs gathers
CCFs_gatherC = zeros(secPT,virTraceNum,virShotNum); 

% the acausal part of all the CCFs gathers
CCFs_gatherAC = zeros(secPT,virTraceNum,virShotNum);
%% Rearrange the original DAS data
% Rearrange the original DAS data according to the CCF duration and overlap
% ratio of adjacent time-window

% the DAS data after rearranging 
rearrData = zeros(secPT,secNum, N); 

for k=1:N
    for j=1:secNum
        beginPt = round((j-1-(j-1)*alp)*secT/dt)+1;
        rearrData(:,j,k) = data(beginPt:beginPt+secPT-1, k);
    end
end
%% cross-correlation: getting all the CCFs gathers
% cross-correlation stack based on the phase-weighted stack (PWS)
window_T = 1.0*dt;              
dataTempC = zeros(secPT,secNum); % causal part of the CCFs, temp variable
dataTempAC = zeros(secPT,secNum);% acausal part of the CCFs, temp variable

% getting all the CCFs gathers by loop
for k=1:virShotNum
    for i=k+1:k+virTraceNum
        for j=1:secNum
            temp = xcorr(rearrData(:,j,i),rearrData(:,j,k));
            dataTempC(:,j) = temp(secPT:end);
            dataTempAC(:,j) = temp(secPT:-1:1);
            
            % normalization
            dataTempC(:,j) = dataTempC(:,j)./max(abs(dataTempC(:,j)));
            dataTempAC(:,j) = dataTempAC(:,j)./max(abs(dataTempAC(:,j)));
        end
        % it means linear stack (LS) about the following two lines
%         dataTempC = sum(dataTempC,2);
%         dataTempAC = sum(dataTempAC,2);
        
        % phase-weighted stack (PWS)
        c = CalcPWSCoefficient(dataTempC,window_T, dt);
        v = 1.0;
        c = c.^(v);
        dataTempC = sum(dataTempC,2).*c;
        c = CalcPWSCoefficient(dataTempAC,window_T, dt);
        v = 1.0;
        c = c.^(v);
        dataTempAC = sum(dataTempAC,2).*c;
        
        % causal and acausal part of the CCFs gather of certain sub-window
        CCFsC(:,i-k) = dataTempC;
        CCFsAC(:,i-k) = dataTempAC;
    end
    CCFs_gatherC(:,:,k) = CCFsC;
    CCFs_gatherAC(:,:,k) = CCFsAC;
end
%% Dispersion Measurement for all the CCFs gathers
r = 1:59;
RangeFreq = [1,40];
ParaVelocity = [50 2000 1];
mId = 'basic-function';
tensorName = 'RR';

k = 1;
[E,freq,v] = CLPhaseShiftofSW(CCFs_gatherC(:,:,k),dt,r,RangeFreq,ParaVelocity,mId,tensorName);

[rowNum,colNum] = size(E);
E_gather = zeros(rowNum,colNum,virShotNum);
E_gather(:,:,k) = E;
for k=2:virShotNum
    [E,freq,v] = CLPhaseShiftofSW(CCFs_gatherC(:,:,k),dt,r,RangeFreq,ParaVelocity,mId,tensorName);
    E_gather(:,:,k) = E;
end

%% plot the dispersion image of the first sub-window
gl = 10;   % gauge-length of the DAS data
figure;
pcolor(freq,v,E_gather(:,:,1));shading interp;
c = colorbar;caxis([0,1]);xlabel(c,'Normalized amplitude');
hold on;
line([ParaVelocity(1)/gl ParaVelocity(2)/gl ],[ParaVelocity(1) ParaVelocity(2)],'Color','black','LineStyle','--');
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');
title('Dispersion image of the first sub-window');