function x_series = trajectory_gen(sys,time_horizon,acc)
% trajectory_gen: function for generating set-points for the case study
%   Code for the Paper entitled "Hierarchical Control for Cyber-Physical Systems via General Approximate Alternating Simulation Relations" in ADHS 2024
% Input:    sys: a gAAS_relation_options object for the system of interest
%           time_horizon: simulation time horizon
%           acc: user-defined parameter for generating the series of set
%           point
% output:   x_series: a series of set-points for x_hat
%   Authors: Bingzhuo Zhong
%   Date: April 1, 2024

    % compute the size of matrices A and B
    x_dim = size(sys.A,2);
    u_dim = size(sys.B,2);
    
    % initialization
    unit = floor(time_horizon/5);
    x_cur = zeros(x_dim,1);
    x_series = zeros(x_dim,time_horizon+1);
    
    % computing the series of desired set-points
    for i = 1:1:time_horizon+1
        x_series(:,i)= x_cur;
        cur_rand = -0.5+rand(1,1);
        if i<= unit
            g = acc+cur_rand;
        elseif unit<i && i<= 2*unit
            g = zeros(u_dim,1)+cur_rand;
        elseif 2*unit<i &&  i<= 4*unit
           g = -acc+cur_rand;
        else
           g = zeros(u_dim,1)+cur_rand;      
        end
        x_nxt = sys.A*x_cur+sys.B*g;
        x_cur = x_nxt;
    end
end