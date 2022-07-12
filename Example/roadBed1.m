%% Modern inversion workflow of the multimodal surface wave DCs
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+
% The code package is expected to solve the mode-losses and mode-aliases of
% the inversion for the multimodal surface wave dispersion curves.
% 
% Idea and operation:
% The pattern search (PS) is used to invert the reliable segment of the 
% fundamental-mode surface wave phase velocities for the first stage. For 
% the second stage, the inverted result of the first stage is set as the 
% initial model, the PS with embedded Kuhn-Munkres (PSEKM) algorithm is 
% adopted for inverting the observed phase velocities of all modes. And for 
% each frequency, a weighted bipartite graph is established between the 
% observed values with no-explicitly-specified-mode-order (NESMO) and 
% predicted values of the model m during the inversion, then the maximum 
% match is determined by the Kuhn-Munkres algorithm for calculating the 
% minimum distance between the observed and predicted data sets. The 
% mode-order information of the observed phase velocities with NESMO would 
% be dynamically evaluated for each model m occurred in the inversion process.
%
% Purpose of this script:
% The script reproduces the inversion result of the roadbed 1 of our
% article.It would generate the Figure 11(c) and 11(d) after runing the
% script.
%
% Note:
% Before using the inversion workflow, we strongly recommend that you read 
% our article and the screen recording tutorial.
% The article and tutorial are placed in "main-folder"\Doc. If you publish
% an article by using our inversion programs, please cite our article. If
% there are any questions about the code, please contact with the first 
% author Yan Yingwei by email "wallace2012y@outlook.com".
%
% Please ensure the main folder is in the current directory of 
% matlab, and its subfolders are also added the path of matlab.
%
% The code has been debugged on the the matlab R2019a. A higher or lower 
% version of matlab may also make the code run successfully.
%
% It would run for about a few minutes.
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
% Revision:  1.0  Date: 5/10/2022
%
% Department of Earth and Space Sciences, Southern University of Science 
% and Technology (SUSTech).
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+

clear;                    % clear the workspace of matlab 
%% load the observed phase velocities
currDir = cd;
filePath = strcat(currDir,'\Data\ObsData.mat');
load(filePath);

% arranging the dispersion-curve (DC) data to get its minimum capacity
pv = arrangeDCData(pv);
[M,N] = size(pv);
%% getting the required information about the fundamental-mode
% the first column is the fundamental-mode
[Index,~] = find(pv(:,1)~=0);
startV = pv(Index(1),1);
endV = pv(Index(end),1);

% the wave-length of the fundamental-mode surface wave
fundamentalLambda = pv(Index,1)./(f(Index)');

% the mean wave-length of the fundamental-mode surface wave
fundamentalLambdaMean = sum(fundamentalLambda)/length(Index); 
%% Inversion parameter settings
num_layer = 4; 

maxIter = 500;
lambda = 1.0;
theta = 0.5;
fTol = 1e-5;
hTol = 0.01;
vTol = 0.1;
misfitType = 'L2';
dStep = [5*ones(1,num_layer) 0.5*ones(1,num_layer-1)];

InvParameterInfo.maxIter = maxIter;
InvParameterInfo.lambda = lambda;
InvParameterInfo.theta = theta;
InvParameterInfo.fTol = fTol;
InvParameterInfo.hTol = hTol;
InvParameterInfo.vTol = vTol;
InvParameterInfo.misfitType = misfitType;
InvParameterInfo.dStep = dStep;
%% Inversion
fValueSequence = zeros(2,maxIter+1);
mInvSequence = zeros(2,2*num_layer-1);
t0Series = zeros(2,1);

% inversion for the first stage, only the fundamental-mode is inverted
tic;
i = 1;
vs_ini = (startV+endV)/0.88/2*ones(1,num_layer);   
h = 1.0*fundamentalLambdaMean/(num_layer-1)*ones(1,num_layer-1);
dBound = [50;50;50*ones(num_layer-2,1);0.5*ones(num_layer-1,1)];
uBound = [1000;1000;1000*ones(num_layer-2,1);10*ones(num_layer-1,1)];

ModelInfo.Bound = [dBound uBound]; % Bound constraints for inversion parameters
ModelInfo.Ini = [vs_ini,h];        % initial model
ModelInfo.vpdvs = 2.45*ones(1,num_layer); % the ratio between vp and vs
ModelInfo.den = 2000*ones(1,num_layer);   % density vector 

ObsInfo.f = f(Index);    % frequency vector of the fundamental-mode
ObsInfo.pv = pv(Index,1);% the observed phase velocity of the fundamental-mode
ObsInfo.maxModeNum = 1;  % the possible minimum mode-order

global n_mode;
n_mode = ObsInfo.maxModeNum;
[mInv,fValue,~] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo);
l = length(fValue);
fValueSequence(i,1:l) = fValue;
mInvSequence(i,:) = mInv;
t0 = toc;
t0Series(i) = t0;

% inversion for the second stage, all the observed phase velocities are
% inverted, the inverted result of the first stage is set as the initial
% model.
i = 2;
ObsInfo.f = f;             % frequency vector
ObsInfo.pv = pv;           % all the observed phase velocities
[~, nn] = size(pv);
ObsInfo.maxModeNum = nn;   % the possible minimum mode-order
n_mode = 12;
ModelInfo.Ini = mInv;
[mInv,fValue,~] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo);
l = length(fValue);
fValueSequence(i,1:l) = fValue;
mInvSequence(i,:) = mInv;
t0 = toc;
t0Series(i) = t0;
%% plot figures
% plot figure 11 (c) in the article
Ind = find(fValueSequence(1,:)~=0);
l1 = length(Ind);
fValue1 = fValueSequence(1,1:l1);
iter1 = 0:l1-1;
Ind = find(fValueSequence(2,:)~=0);
l2 = length(Ind);
fValue2 = fValueSequence(2,1:l2);
iter2 = l1-1:l1+l2-2;
figure;
h1 = plot(iter1,fValue1,'k-',iter2,fValue2,'b-');
axis([0,l1+l2-2 0 90]);
legend(h1,'for the first stage','for the second stage');
xlabel('Iteration');ylabel('RMS (m/s)');

% plot figure 11 (d) in the article
maxDepth = 15;
[vsIniSta, depthIniSta] = calcStairsData(vs_ini,h, maxDepth);
[vsInv1Sta, depthInv1Sta] = calcStairsData(mInvSequence(1,1:num_layer),mInvSequence(1,num_layer+1:end), maxDepth);
[vsInv2Sta, depthInv2Sta] = calcStairsData(mInvSequence(2,1:num_layer),mInvSequence(2,num_layer+1:end), maxDepth);
figure;
h2 = plot(vsIniSta,depthIniSta,'k--',vsInv1Sta,depthInv1Sta,'b-',vsInv2Sta,depthInv2Sta,'r-');
set(gca,'YDir','reverse');
axis([100 900 0 15]);legend(h2,'Initial','for the first stage','for the second stage');
xlabel('S-wave velocity (m/s)');ylabel('Depth (m)');






