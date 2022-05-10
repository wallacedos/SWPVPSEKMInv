function [vsSta, depthSta] = calcStairsData(vs,h, maxDepth)
if nargin==2
   maxDepth = 1.5*sum(h);
end

n = length(vs);
vs = [vs(1) vs];
depth = zeros(1,n+1);
for i=2:n
    depth(i) = sum(h(1:i-1));
end
depth(end) = maxDepth;
[vsSta,depthSta] = stairs(vs,depth);
end

