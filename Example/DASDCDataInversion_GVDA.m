%% Multimodal dispersion inversion of the Long Line I DAS data
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-
% Purpose of this script:
% This is a script for inverting the multimodal dispersion curves of the 
% Long Line I distributed-acoustic-sensing (DAS) data.
% And the DAS data sets are collected at the Garner Valley Downhole 
% Array (GVDA) field test site in Southern California in September 2013, 
% and the survey area is adjacent to Highway 74 of California. 
%
% About the information of survey area and data, please visit:
% https://gdr.openei.org/submissions/481
% http://nees.ucsb.edu/facilities/GVDA
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
% It would run for about ten hours.
%
% Our article:
% Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
% the multimodal surface wave dispersion curves: Staging strategy and Pattern 
% search with embedded Kuhn-Munkres algorithm, Geophysical Journal
% International,231(01), 47-71, 
% https://doi.org/10.1093/gji/ggac178. 
%
% Author(s): Yan Yingwei
% Email:     wallace2012y@outlook.com
% Copyright: 2022-2025 
% Revision:  1.0  Date: 7/11/2022
%
% Department of Earth and Space Sciences, Southern University of Science 
% and Technology (SUSTech).
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-

clear;         % clear the workspace of matlab
%% Inversion settings
ptNum = 125;   % the number of survey points
numLayer = 6;  % the number of layers

maxIter = 300; % maximum iterations
lambda = 1.0;  % expand factor
theta = 0.5;   % constriction factor
fTol = 0.01;   % tolerance for misfit function value

% tolerance for layer thickness, would be adaptively determined in the 
% following part
hTol = 0.2; 

% tolerance for layer S-wave velocity, would be adaptively determined in the 
% following part
vTol = 1;   

% the type of the misfit function, 'L2' means deriving by the L2-norm
misfitType = 'L2';

% initial step-length for updating the model parameter (S-wave velocity and
% layer thickness). It would be adaptively determined in the following part
dStep = [10*ones(1,numLayer) 2*ones(1,numLayer-1)];

InvParameterInfo.maxIter = maxIter;
InvParameterInfo.lambda = lambda;
InvParameterInfo.theta = theta;
InvParameterInfo.fTol = fTol;
InvParameterInfo.hTol = hTol;
InvParameterInfo.vTol = vTol;
InvParameterInfo.misfitType = misfitType;
InvParameterInfo.dStep = dStep;

%% Inversion for obtaing the 1D S-wave velocity of all sub-windows 
% inverted 1D S-wave velocity of all sub-windows 
mInvSequence = zeros(2*ptNum,2*numLayer-1);

% the misfit function values of the whole inversion
fValueSequence = zeros(2*ptNum,maxIter+1);

% computing time information of the whole inversion
t0 = zeros(2*ptNum,1);

currDir = cd;
filePath = strcat(currDir,'\Data\');
for i=1:ptNum
    curFileName = strcat(filePath,'ObsData_GVDA_',num2str(i),'.mat');
    load(curFileName);
    [M,N] = size(pv);
    
    [Index,~] = find(pv(:,1)~=0);
    startV = pv(Index(1),1);
    endV = pv(Index(end),1);
    
    % wavelength vector of fundamental surface-wave
    fundamentalLambda = pv(Index,1)./(f(Index)');
    
    % maximum wavelength of surface waves
    fundamentalLambdaMax = fundamentalLambda(1);
    
    % constructing the initial model according to the fundamental-mode
    % surface wave dispersion curve
    % initial S-wave velocity
    vsIni = (startV+endV)/0.88/2*ones(1,numLayer);
    % initial layer thickness vector
    hIni = 0.5*fundamentalLambdaMax/(numLayer-1)*ones(1,numLayer-1);
    
    % the following parameters are adaptively determined
    InvParameterInfo.dStep = [vsIni(1)/20*ones(1,numLayer) hIni(1)/20*ones(1,numLayer-1)];
    InvParameterInfo.hTol = hIni(1)/500;
    InvParameterInfo.vTol = vsIni(1)/500;
    
    vsUBound = max(max(pv))/0.88*2.5;
    vsDBound = 80;
    Bound = [vsDBound*ones(numLayer,1),vsUBound*ones(numLayer,1);0.005*fundamentalLambdaMax*ones(numLayer-1,1),1.5*fundamentalLambdaMax*ones(numLayer-1,1)];
    
    % the ratio between the P-wave velocity and S-wave velocity, invariant
    % during the inversion process
    vpdvs = 3.0*ones(1,numLayer);
    
    % first stage inversion: only inverting the fundamental-mode surface
    % wave dispersion curve
    ModelInfo.Bound = Bound;
    ModelInfo.Ini = [vsIni,hIni];
    ModelInfo.vpdvs = vpdvs;
    
    % density of the model, invariant during the inversion process
    ModelInfo.den = 2000*ones(1,numLayer);
    ObsInfo.f = f(Index);     % frequency vector of the fundamental-mode
    ObsInfo.pv = pv(Index,1); % observed values of the fundamental-mode
    ObsInfo.maxModeNum = 1;   % the possible minimum mode-order!!!
    
    global n_mode;
    n_mode = ObsInfo.maxModeNum;
    fprintf('The %d-th survey point, first stage...\n',i);
    tic;
    [mInv,fValue, ~] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo);
    t0(2*i-1) = toc;
    mInvSequence(2*i-1,:) = mInv;
    l = length(fValue);
    fValueSequence(2*i-1,1:l) = fValue;
    
    % second stage inversion: inverting all the phase velocities
    % the initial model is set as the inversion result of the first stage
    [~,colNum] = size(pv);
    ModelInfo.Ini = mInv;
    ObsInfo.f = f;
    ObsInfo.pv = pv;
    ObsInfo.maxModeNum = colNum;
    
    n_mode = 50;
    fprintf('The %d-th survey point, second stage...\n',i);
    tic;
    [mInv,fValue, ~] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo);
    t0(2*i) = toc;
    mInvSequence(2*i,:) = mInv;
    l = length(fValue);
    fValueSequence(2*i,1:l) = fValue;
end

% the final misfit value of the first stage inversion
fValue1End = zeros(1,ptNum);  
% the final misfit value of the second stage inversion
fValue2End = zeros(1,ptNum);

for i=1:ptNum
    Ind = find(fValueSequence(2*i-1,:)~=0);
    fValue1End(i) = fValueSequence(2*i-1,Ind(end));
    Ind = find(fValueSequence(2*i,:)~=0);
    fValue2End(i) = fValueSequence(2*i,Ind(end));
end
%% Kriging interpolation on inversion results at the second stage
% it would generate the pseudo-2D S-wave velocity profile
Index = find(fValue2End<12);
lIndex = length(Index);

X = zeros(lIndex*numLayer,1);  % horizontal position (m)
Y = zeros(lIndex*numLayer,1);  % vertical position (m), depth 
Z = zeros(lIndex*numLayer,1);  % inverted S-wave velocity values

for i=1:lIndex
    for j=1:numLayer
        X((i-1)*numLayer+j) = Index(i)-1;
        Y((i-1)*numLayer+j) = sum(mInvSequence(2*Index(i),numLayer+1:numLayer+j-1));
        Z((i-1)*numLayer+j) = mInvSequence(2*Index(i),j);
    end
end

[VSTomo,gridX,gridY] = kriging(X,Y,Z);
VSTomo = VSTomo';
sigma = 20;
VSSmoothTomo = imgaussfilt(VSTomo,sigma); % gaussian smooth filter
%% plot the pseudo-2D S-wave velocity profile of the survey area
% the pseudo-2D S-wave velocity profile generated from the second-stage
% inverion results
figure;
pcolor(gridX+231,gridY, VSSmoothTomo);shading interp; axis image;
set(gca,'YDir','reverse');
colormap(hsv);caxis([100 1200]);c = colorbar;xlabel(c,'S-wave velocity (m/s)');
set(c,'limits',[100 1000]); ylim([0, 120]);
xlabel('Position (m)');ylabel('Depth (m)');
title('The pseudo-2D S-wave velocity profile of the second-stage');

