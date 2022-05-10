function match = KM(r)
% Summary of this function goes here.
% match = KM(r)
% Detailed explanation goes here.
% The function is for getting the maximum matching of weighted bipartite
% graphs. 
%   IN      
%      r: the matrix used to describe the bipartite graph.             
%
%  OUT   
%      match: the output matching relationship for the bipartite graph.
%
%  EXAMPLE
%      r = [3,0,10,2;2,1,3,2;0,0,5,9];
%     match = KM(r);
% 
%     match would be [3;1;4], it means that the 3rd element of the first
%     row, the 1st element of the second row,the 4th element of the third
%     row are selected, is the maximum match.
%
%  References: 
%  Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
%  the multimodal surface wave dispersion curves: Staging strategy and Pattern 
%  search with embedded Kuhn-Munkres algorithm, Geophysical Journal
%  International (online),
%  https://doi.org/10.1093/gji/ggac178. 
%
%  Kuhn, H. 1955. The Hungarian method for the assignment problem, Naval 
%  Research Logistics, 2, 83-97, https://doi.org/10.1002/nav.3800020109.
%   
%  West, D. 2020. Introduction to Graph Theory Second Edition (in Chinese), 
%  China Machine Press.
% 
%  Author(s): Yan Yingwei
%  Copyright: 2022-2025 
%  Revision:  1.0  Date: 5/10/2022
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

global g_map;
global g_m;
global g_n;
global g_label_value_x;
global g_label_value_y;
global g_match_x;
global g_match_y;
global g_min_weight_diff;
global g_visible_x;
global g_visible_y;
%%
%init
[g_m,g_n]=size(r);
g_map = r;
g_label_value_x = zeros(g_m,1);
g_label_value_y = zeros(g_n,1);
g_match_x = ones(g_m,1)*-1;
g_match_y = ones(g_n,1)*-1;
for i=1:g_m
    for j=1:g_n
        if g_label_value_x(i)<g_map(i,j)
            % Assign the vertex label to the vertices of set X
            g_label_value_x(i) = g_map(i,j);
        end
    end
end
%%
%KM-algorithm
for i =1:g_m
    while(1)
        g_min_weight_diff=10^10;
        g_visible_x = zeros(g_m,1);
        g_visible_y = zeros(g_n,1);
        if dfs(i)
            break;
        end
        for j=1:g_m
            if g_visible_x(j)
                g_label_value_x(j)=g_label_value_x(j)-g_min_weight_diff;
            end
        end
        for j=1:g_n
            if g_visible_y(j)
                g_label_value_y(j)=g_label_value_y(j)+g_min_weight_diff;
            end
        end
    end
end
match = g_match_x;
end

