function [mInv,fValueSequence,mAccepted] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo)
% Summary of this function goes here.
% [mInv,fValueSequence, mAccepted] = rayleighDCMNKPSInv(ObsInfo,ModelInfo,InvParameterInfo)
% Detailed explanation goes here.
% The function is for inverting the fundamental/multimodal Rayleigh wave 
% dispersion curves by the pattern search with embeded Kuhn-Munkres (PSEKM) 
% algorithm.
%
%   IN      
%        ObsInfo: the struct about the observed information, it consists of
%                 following fields, 'f' (frequency axis vector), 'pv'
%                 (observed phase velocity values matrix, the first column 
%                 represents fundamental-mode observations, the mode-order
%                 of the other observed phase velocities would be determined,
%                 adaptively, 0-elements are the unobserved data point),
%                'maxModeNum', denotes the possible max mode-order. 
%      ModelInfo: the struct about the model information during the inversion,
%                 it consists of following fields,'Bound' (searching bound
%                 for each paramter, it's n*2 matrix, n is the number of
%                 reconstructed parameter), 'Ini' (initial model of the inverion)
%                 'den' (density of the inversion, it remains a constant),
%                 'vpdvs' (the ratio of P-wave and S-wave velocity for each
%                 layer, it's also a row vector.).
%InvParameterInfo:it's a struct about the inversion parameter setting, it 
%                 consists of following fields, 'dStep' (initial step-length
%                 vector for all the inversion parameters), 'lambda', (the
%                 expansion factor for the pattern search, generally greater 
%                 than or equal to 1), 'theta',(constraction factor,
%                 generally less than 1), 'fTol', (the tolerance parameter
%                 of the misfit function, when the misfit value is less than
%                 fTol times the initial misfit, the iteration is
%                 terminated), 'hTol',(tolerance parameter for the
%                 layer thickness), 'vTol',(tolerance parameter for the
%                 S-wave velocity),  'maxIter' (the maximum iterations of 
%                 inversion),'misfitType' (the choice of type for misfit 
%                 function, it consists of 'L1', 'L2', 'L2' is recommended.             
%
%  OUT   
%          mInv:  the inverted result after iterations.
% fValueSequence: the misfit function values of the iterations.
%      mAccepted: the accepted models of each iteration.
%
%  EXAMPLE
%  Please read the roadBed1.m at the path "Example".
%
%
%  References: 
%  Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
%  the multimodal surface wave dispersion curves: Staging strategy and Pattern 
%  search with embedded Kuhn-Munkres algorithm, Geophysical Journal
%  International (online),
%  https://doi.org/10.1093/gji/ggac178. 
%   
%
%  Author(s): Yan Yingwei
%  Copyright: 2022-2025 
%  Revision:  1.0  Date: 5/10/2022
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

%% read the input parameters
freq = ObsInfo.f;                     % frequency-vector (Hz), a row vector
pvObserved = ObsInfo.pv;              % the observed values, a matrix
maxModeNum = ObsInfo.maxModeNum;      % the possible maximum mode-order

ModelBoundConstraint = ModelInfo.Bound; 
mCurrIni = ModelInfo.Ini; % initial model, a row vector, it's like [vs, h]
vpdvs = ModelInfo.vpdvs;  % the ratio of vp and vs for each layer
den = ModelInfo.den; % density vector, it's a constant during the inversion

dStep = InvParameterInfo.dStep; 
maxIter = InvParameterInfo.maxIter;
lambda = InvParameterInfo.lambda;
theta =InvParameterInfo.theta;
fTol = InvParameterInfo.fTol;
hTol = InvParameterInfo.hTol;
vTol = InvParameterInfo.vTol;
misfitType = InvParameterInfo.misfitType; 
%% inversion by iterations
parameterNum = length(mCurrIni); % the number of inverted parameters
n = length(den);                 % the number of layer

B=eye(parameterNum,parameterNum);
p=2*parameterNum;
C=zeros(parameterNum,parameterNum);
C(:,1:parameterNum)=B(:,1:parameterNum);
C(:,parameterNum+1:p)=-B(:,1:parameterNum);
D=B*C;

fValueSequence = zeros(maxIter+1,1);
fCurrSearchValue = zeros(p,1);
pvPredicted = calcmulti(freq,mCurrIni(1:n),mCurrIni(n+1:end), vpdvs.*mCurrIni(1:n),den);

fCurrValue = calcObjfOfKM(pvObserved,pvPredicted,misfitType);
fValueSequence(1) = fCurrValue;

% the expected misfit value that can be achieved after the inversion.
fExpectedValue = fTol*fCurrValue;    

% the current search models of each iteration
mCurrSearch = zeros(p,parameterNum); 

% the accepted models for the whole inversion
mAccepted = zeros(maxIter,parameterNum);

% iteration
ii = 0;
for i=1:maxIter
    fCurrSearchValue(:) = inf;
    mCurrSearch(:,:) = 0;
    for j=1:p
        ii = ii+1;
        kk = ceil(mod(j,parameterNum+0.1));
        mNextSearch = mCurrIni+dStep(kk)*D(:,j)'; % search model
        Bounds =  ModelBoundConstraint(kk,:);
        
        if(mNextSearch(kk)<Bounds(1))
            mNextSearch(kk) = Bounds(1);         % rebound
        end
        
        if(mNextSearch(kk)>Bounds(2))
            mNextSearch(kk) = Bounds(2);        % rebound
        end
        
        try
            pvPredicted = calcmulti(freq,mNextSearch(1:n),mNextSearch(n+1:end), vpdvs.*mNextSearch(1:n),den);
            [~,nn] = size(pvPredicted);
            if(nn>=maxModeNum)
                fCurrSearchValue (j) = calcObjfOfKM(pvObserved,pvPredicted,misfitType);
                mCurrSearch(j,:) = mNextSearch;
            end
        catch
            continue;
        end
    end
    
    [~,I] = min(fCurrSearchValue);
    
    if fCurrValue/fCurrSearchValue(I)>1
        fCurrValue = fCurrSearchValue(I);
        mCurrIni =  mCurrSearch(I,:);
        dStep = lambda*dStep;
    else
        dStep = theta*dStep;
    end
    
    mAccepted(i,:) = mCurrIni;
    
    fValueSequence(i+1) = fCurrValue;
    fprintf('The %d iteration current objective function value is %f.\n',i,fCurrValue);
    
    % the iteration is terminated when the following conditions are
    % satisfied.
    if(fCurrValue<fExpectedValue || min(dStep(1:n))<vTol || min(dStep(n+1:end))<hTol)       
        break;
    end
end

fValueSequence = fValueSequence(1:i+1);
[~,I] = min(fValueSequence);
mAccepted = mAccepted(1:i,:);
if(I==1)
    mInv = mAccepted(I,:);
else
    mInv = mAccepted(I-1,:);
end
end

