% check_gamma1_R_script: MATLAB script for computing gamma1 and matrix R using the results in (Girard and Pappas, 2009, Proposition 1), see equation (16)
% Code for the Paper entitled "Hierarchical Control for Cyber-Physical Systems via General Approximate Alternating Simulation Relations" in ADHS 2024
%   Authors: Bingzhuo Zhong
%   Date: April 1, 2024

% ATTENTION: Please run ADHS_REP_gAAS.m before running the current script

clear all

% loading parameters for the case study of the vehicle
load('ADHS_2D_Result.mat');

% input bound for the abstract system
uhat_bd = 35;

% The following results are reported in the case study of the paper

% compute R accroding to (Girard and Pappas, 2009, Proposition 1), see equation (16) in the paper
R = (car.B'*car.M*car.B)^(-1)*car.B'*car.M*car.P*car.B_hat

% compute gamma1_pri accroding to (Girard and Pappas, 2009, Proposition 1), see equation (16) in the paper
gamma1_pri = sqrt((car.B*R-car.P*car.B_hat)'*car.M*(car.B*R-car.P*car.B_hat)*uhat_bd*uhat_bd)