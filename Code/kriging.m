function [elevation,gridX,gridY] = kriging(x,y,z,range,sill)
% KRIGING  Ordinary Kriging Interpolation
%
%   3D interpolation of scattered height data against x and y data. 
%
%   USES:
%   [elevation]                 = KRIGING(PointsX,PointsY,PointsElev)
%   [elevation]                 = KRIGING(PointsX,PointsY,PointsElev,Range,Sill)
%   [elevation,gridX,gridY]     = KRIGING(PointsX,PointsY,PointsElev)
%   [elevation,gridX,gridY]     = KRIGING(PointsX,PointsY,PointsElev,Range,Sill)
%   
%   The input variables range and sill are optional; 'Default' values will be
%   used if not supplied. The 'default' values are not correct but give a
%   good indiciation of your kriging. To determine the right values use
%   variography.
% 
%   See also SURF, MESHGRID.
%
%   Author: J.W. Buist
%   Data:   15-5-2016

%% Checking input
if nargin < 3
    error('Error. Input at least PointsX, PointsY and PointsElev')
end

if ~exist('range','var')
    range = 26440.092;
end

if ~exist('sill','var')
    sill = 62583.893;
end

%% Calculating trend
xx = x(:);
yy = y(:);
zz = z(:);
N = length(xx);
O = ones(N,1);
C = [xx yy O]\zz;
PointsX = x;
PointsY = y;
PointsElev = z - (C(1).*x + C(2).*y +C(3));
L = length(PointsX);
S = size(PointsElev);
if S(1) > S(2)
    PointsElev = PointsElev.';
end

%% Kriging Interpolation
% Building grid
gridX = linspace(min(PointsX),max(PointsX),L);
gridY = linspace(min(PointsY),max(PointsY),L);

% Building reduncy matrix
ReduMatrix = zeros(L,L);
for i = 1:L
    for j = 1:L
        Distance = sqrt((PointsX(i) - PointsX(j))^2 + (PointsY(i) - PointsY(j))^2);
        if Distance > range
            ReduMatrix(i, j) = sill;
        else
            ReduMatrix(i, j) = sill*((3*Distance./(2*range)) -1/2*(Distance ./range).^3);
        end
    end
end

% Kriging Interpolation
xCords = zeros(L,L);
yCords = zeros(L,L);
elevation = zeros(L,L);
hw = waitbar(0,'Kriging...','CreateCancelBtn', {@(H,~) delete( findobj( get(H,'Parent'), '-Depth', 1, 'Tag', 'TMWWaitbar'))});
for gridXcor = 1:L
    try
        waitbar(gridXcor/L,hw)
    catch
        error('Kriging cancelled by user')
    end
    
    for gridYcor = 1:L
        ProxVector = zeros(L, 1);
        for a = 1:L
            Distance = sqrt((PointsX(a) - gridX(gridXcor))^2 + (PointsY(a) - gridY(gridYcor))^2);
           if Distance > range
            ProxVector(a) = sill;
           else
            ProxVector(a) = sill*((3*Distance./(2*range))-1/2*(Distance ./range).^3);
           end
        end
        
        Weights = ReduMatrix \ ProxVector;
        XYElev = PointsElev * Weights;
        XYElev = XYElev + ((C(1) * gridX(gridXcor)) + (C(2) * gridY(gridYcor)) + C(3) );
        xCords(gridXcor,gridYcor) = gridX(gridXcor);
        yCords(gridXcor,gridYcor) = gridY(gridYcor);
        elevation(gridXcor,gridYcor) = XYElev;    
    end
end
delete(hw)