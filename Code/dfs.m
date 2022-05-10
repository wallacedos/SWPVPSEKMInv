function res=dfs(i)
% Summary of this function goes here.
% res=dfs(i)
% Detailed explanation goes here.
% The function implements the depth-first search (DFS).
% It is used in the process of using the Kuhn-Munkres algorithm to solve
% the maximum matching of the weighted bipartite graph.

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

g_visible_x(i) = 1;
for j =1:g_n
    if(~g_visible_y(j))
        tmp = g_label_value_x(i)+g_label_value_y(j)-g_map(i,j);
        if tmp<10^-5
            g_visible_y(j)=1;
            if(g_match_y(j)==-1||dfs(g_match_y(j)))
                g_match_x(i)=j;
                g_match_y(j)=i;
                res=true;
                return
            end
        elseif tmp>0
            if tmp<g_min_weight_diff
                g_min_weight_diff  =tmp;
            end
        end
    end
end
res = false;
end
