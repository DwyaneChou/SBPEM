% Choose the spatial discrete scheme for split
% pass = 0 for the full L_operator
% pass = 1 for the fast pass operator
% pass = 2 for the slow pass operator
function [LU,LV,LZ] = spatial_discrete(pass,STATE,MESH)
                             
if pass == 0
    [LU,LV,LZ] = L_operator(STATE,MESH);
elseif pass == 1
    [LU,LV,LZ] = fast_pass(STATE,MESH);
elseif pass == 2
    [LU,LV,LZ] = slow_pass(STATE,MESH);
end