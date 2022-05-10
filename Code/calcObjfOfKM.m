function fValue = calcObjfOfKM(pvObserved,pvPredicted, misfitType)
% Summary of this function goes here.
% fValue = calcObjfOfKM(pvObserved,pvPredicted, misfitType)
% Detailed explanation goes here.
% The function is for measuring the misfit value between the observed and 
% predicted phase velocities.The mode-order information would be evaluated,
% adaptively. The mode-oder is got by the Kuhn-Munkres algorithm to solve
% the bipartite graph. 
%   IN      
%     pvObserved: the observed phase velocities, it's a matrix, the first
%                 column represents the fundamental-mode, the mode-order of
%                 the other observed values would be evaluated, adaptively.
%    pvPredicted: the predicted phase velocities.
%     misfitType: the choice of type for misfit function, it consists of 'L1', 
%                 'L2', 'L2' is recommended.             
%
%  OUT   
%         fValue: the misfit function value between the observed and
%                 predicted phase velocities.
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

if nargin==2
    misfitType = 'L2';
end
fValue = 0;
[M,N] = size(pvObserved);

if strcmp(misfitType,'L2')
    [Ind,~] = find(pvObserved(:,1)~=0);
    if N==1  % only the fundamental-mode is inverted
        temp = (pvObserved(Ind,1)-pvPredicted(Ind,1));
        temp = temp.*temp;
        fValue = fValue+sum(temp);
        dataPointNum = length(Ind);
        fValue = sqrt(fValue/dataPointNum);
    else
        dataPointNum = 0;
        temp = (pvObserved(Ind,1)-pvPredicted(Ind,1));
        temp = temp.*temp;
        fValue = fValue+sum(temp);
        dataPointNum = dataPointNum+length(Ind);
        for i=1:Ind(1)-1
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -(pvObserved(i,j)-pvPredicted(i,1:end)).*(pvObserved(i,j)-pvPredicted(i,1:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue = fValue+(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)))*(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)));
            end
            clear r Index;
        end
        for i=Ind(1):Ind(end)
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -(pvObserved(i,j)-pvPredicted(i,2:end)).*(pvObserved(i,j)-pvPredicted(i,2:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue = fValue+(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)+1))*(pvObserved(i, Index(ii))-pvPredicted(i,match(ii)+1));
            end
            clear r Index;
        end
        for i=Ind(end)+1:M
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -(pvObserved(i,j)-pvPredicted(i,1:end)).*(pvObserved(i,j)-pvPredicted(i,1:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue =fValue+(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)))*(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)));
            end
            clear r Index;
        end
        fValue = sqrt(fValue/dataPointNum);
    end
elseif strcmp(misfitType,'L1')
    [Ind,~] = find(pvObserved(:,1)~=0);
    if N==1        % only the fundamental-mode is inverted
        temp = abs(pvObserved(Ind,1)-pvPredicted(Ind,1));
        fValue = fValue+sum(temp);
        dataPointNum = length(Ind);
        fValue = fValue/dataPointNum;
    else
        dataPointNum = 0;
        temp = abs(pvObserved(Ind,1)-pvPredicted(Ind,1));
        fValue = fValue+sum(temp);
        dataPointNum = dataPointNum+length(Ind);
        for i=1:Ind(1)-1
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -abs(pvObserved(i,j)-pvPredicted(i,1:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue = fValue+abs(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)));
            end
            clear r Index;
        end
        for i=Ind(1):Ind(end)
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -abs(pvObserved(i,j)-pvPredicted(i,2:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue = fValue+abs(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)+1));
            end
            clear r Index;
        end
        for i=Ind(end)+1:M
            k = 0;
            for j=2:N
                if(pvObserved(i,j)~=0)
                    k = k+1;
                    Index(k) = j;
                    r(k,:) = -abs(pvObserved(i,j)-pvPredicted(i,1:end));
                    dataPointNum = dataPointNum+1;
                end
            end
            if exist('r','var')
                match = KM(r);
            end
            for ii=1:k
                fValue = fValue+abs(pvObserved(i,Index(ii))-pvPredicted(i,match(ii)));
            end
            clear r Index;
        end
        fValue = fValue/dataPointNum;
    end
end
end


