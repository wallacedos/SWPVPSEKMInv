function pvObsA = arrangeDCData(pvObs)
% Summary of this function goes here.
% pvObsA = arrangeDCData(pvObs)
% Detailed explanation goes here.
% The function is for arranging the disperion-curve (DC) data to get its
% minimum capacity.
%
%  Author(s): Yan Yingwei
%  Copyright: 2022-2025 
%  Revision:  1.0  Date: 5/10/2022
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

[~,N] = size(pvObs);
for j=3:N
    [Index,~] = find(pvObs(:,j)~=0);
    for k=Index'
        for ii=2:j-1
            if(pvObs(k,ii)==0)
                pvObs(k,ii) = pvObs(k,j);
                pvObs(k,j) = 0;
                break;
            end
        end
    end
end
for j=N:-1:1
    [Index,~] = find(pvObs(:,j)~=0);
    if ~isempty(Index)
        break;
    end
end
pvObsA = pvObs(:,1:j);
end

