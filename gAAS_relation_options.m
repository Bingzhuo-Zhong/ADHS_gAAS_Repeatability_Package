classdef gAAS_relation_options
    %   gAAS_relation_options: options for computing gAAS relation 
    %   Code for the Paper entitled "Hierarchical Control for Cyber-Physical Systems via General Approximate Alternating Simulation Relations" in ADHS 2024
    %   Authors: Bingzhuo Zhong
    %   Date: April 1, 2024
    
    properties
        % system dynamics for the concrete and abstract system
        A;              % A matrix for the concrete system
        B;              % B matrix for the concrete system
        C;              % C matrix for the concrete system
        A_hat;          % A matrix for the abstract system
        B_hat;          % B matrix for the abstract system
        C_hat;          % C matrix for the abstract system
        X_hat;          % State space of the abstract system (a cell containing multiple polytopes)
        U_hat;          % Input space of the abstract system (a cell containing multiple polytopes)
        % marices for the model order reduction
        P;              % P matrix for model order reduction
        D;              % D matrix for the geometric condition
        Q;              % Q matrix for the geometric condition
        % matrices for the general Approximate Alternating Simulation Relation
        K;              % K matrix for the interface function
        S;              % S matrix for the interface function
        R;              % R matrix for the interface function
        M;              % M matrix for the simulation relation
        varepsilon;     % Constant epsilon for the simulation relation
        % parameters for computing the simulation relation
        kappa;          % Parameter kappa in equation (13)
        M_bar;          % Matrix M_bar in equation (13) and (14)
        K_bar;          % Matrix K_bar in equation (13) and (14)
        gamma1;         % gamma1 in equation (11)
        gamma2;         % gamma2 in equation (12)
        opt_option;     % option for solving the simulation relation
    end
    
    methods
        function obj = gAAS_relation_options()
            % gAAS_relation_options: Instantiating an object for computing gAAS relations
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % initialization 
            obj.opt_option.solver = 'mosek'; % installation of Mosek required. Any SDP supported by YALMIP can be deployed, see https://yalmip.github.io/allsolvers/
            obj.opt_option.verbose = 0;      % not printing the details for solving the optimization problem.
        end
        
        %% Computation of the gAAS relation: compute K and M jointly
        function [fea,obj] = solve_kappa(obj,kappa)
            % solve_kappa: solving M and K jointly according to the conditions in (13) and (14)
            % Input:    obj: a gAAS_relation_options object
            %           kappa: parameter kappa in equation (13)
            % output:   fea: flag indicating whether the optimization problem is feasible 
            %               fea = 1: the optimization problem is feasible
            %               fea = 0: the optimization problem is not feasible
            %           obj: a gAAS_relation_options object with the computed simulation relation (if feasible) 
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % computing the sizes of matrices A, B, and C
            size_x = size(obj.A,2);
            size_u = size(obj.B,2);
            size_y = size(obj.C,1);
            
            % define variables to be solved
            M_b = sdpvar(size_x,size_x);
            K_b = sdpvar(size_u,size_x);
            
            % specifying constraints as in equations (13) and (14)
            con1 = [kappa*M_b M_b'*obj.A'+K_b'*obj.B';obj.A*M_b+obj.B*K_b  M_b]>=0;
            con2 = [M_b M_b*obj.C';obj.C*M_b eye(size_y)];
            con3 = M_b>=0;
            constraints = [con1,con2,con3];
            
            % defining the objective function
            fun = -logdet(M_b);
            
            % solving the SDP problem
            ops = sdpsettings('solver', obj.opt_option.solver, 'verbose',obj.opt_option.verbose);
            sol = optimize(constraints, fun, ops);
            
            if sol.problem == 0
                % solve sucessfully
                
                % saving the current results
                obj.kappa = kappa;
                obj.M = inv(value(M_b));
                obj.K =  value(K_b)*obj.M;
                obj.M_bar = value(M_b);
                obj.K_bar = value(K_b);
                
                % compute matrices S and R in the interface function in (10)
                % according to equation (11) and (17)
                [obj,fea2] = obj.compute_interface();
                
                % check feasibility condition in equation (15)
                if fea2 == 1
                    [obj,fea3] = obj.check_fea();
                else
                    fea3 = 0;
                end
                
                % report the final results 
                if fea2 == 1 && fea3 == 1
                    fea = 1;
                    disp('A general Approximate Alternating Simulation Relation can be found.');
                else
                    fea = 0;
                    disp('General Approximate Alternating Simulation Relation can not be found.');
                end
            else
                % the optimization problem is not feasible
                fea = 0;
                obj.kappa = kappa;
                obj.M = [];
                obj.K = [];
                obj.M_bar = [];
                obj.K_bar = [];
                obj.R = [];
                obj.S = [];
                disp('General Approximate Alternating Simulation Relation can not be found.');
            end
            
        end

       %% Computation of the gAAS relation: fix K and compute M
        function [fea,obj] = solve_kappa_fixK(obj,kappa,K)
            % solve_kappa_fixK: solving M according to the condition in (13) and (14) by fixing K
            % Input:    obj: a gAAS_relation_options object
            %           kappa: parameter kappa in equation (13)
            %           K: matrix K in the interface function in equation (10)
            % output:   fea: flag indicating whether the optimization problem is feasible 
            %               fea = 1: the optimization problem is feasible
            %               fea = 0: the optimization problem is not feasible
            %           obj: a gAAS_relation_options object with the computed simulation relation (if feasible) 
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            %  computing the sizes of matrices A and C
            size_x = size(obj.A,2);
            size_y = size(obj.C,1);
            
            % define variables to be solved
            M_b = sdpvar(size_x,size_x);
            K_b = K*M_b;
            
            % specifying constraints as in equations (13) and (14)
            con1 = [kappa*M_b M_b'*obj.A'+K_b'*obj.B';obj.A*M_b+obj.B*K_b  M_b]>=0;
            con2 = [M_b M_b*obj.C';obj.C*M_b eye(size_y)];
            con3 = M_b>=0;
            constraints = [con1,con2,con3];
            
             % defining the objective function
            fun = -logdet(M_b);
            
            % solving the SDP problem
            ops = sdpsettings('solver', obj.opt_option.solver, 'verbose',obj.opt_option.verbose);
            sol = optimize(constraints, fun, ops);
            
            if sol.problem == 0
                % solve sucessfully
                % saving the current results
                obj.kappa = kappa;
                obj.M = inv(value(M_b));
                obj.K = K;
                obj.M_bar = value(M_b);
                obj.K_bar = [];
                
                % compute matrix S and R in the interface function in (10)
                % according to equation (11) and (17)
                [obj,fea2] = obj.compute_interface();
                
                % check feasibility condition in equation (15)
                if fea2 == 1
                    [obj,fea3] = obj.check_fea();
                else
                    fea3 = 0;
                end
                
                % report the final results 
                if fea2 == 1 && fea3 == 1
                    fea = 1;
                    disp('A general Approximate Alternating Simulation Relation can be found.');
                else
                    fea = 0;
                    disp('General Approximate Alternating Simulation Relation can not be found.');
                end
            else
                % the optimization problem is not feasible
                fea = 0;
                obj.kappa = kappa;
                obj.M = [];
                obj.K = K;
                obj.M_bar = [];
                obj.K_bar = [];
                obj.R = [];
                obj.S = [];
                disp('General Approximate Alternating Simulation Relation can not be found.');
            end
        end
        
        %% Computation of the gAAS relation: computing S and R in the interface function give K and M
        function [obj,fea] = compute_interface(obj)
            % compute_interface: computing the interface function given marices K and M
            % Input:    obj: a gAAS_relation_options object
            % output:   fea: flag indicating whether matrices S and R can be computed
            %               fea = 1: feasible S and R in founded
            %               fea = 0: feasible S or/and R in not founded due to unknown computational issue
            %           obj: a gAAS_relation_options object with the computed interface function (if feasible) 
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % compute dimensions of the state, input, and output space of the concrete and the abstract system
            size_x = size(obj.A,2);
            size_u = size(obj.B,2);
            size_y = size(obj.C,1);
            size_uhat = size(obj.B_hat,2);
            
            % initializing variables for computing matrices S and R
            S_c = sdpvar(size_x,size_uhat);
            R_c = sdpvar(size_u,size_uhat);
            
            % formulate equation (17) as a constraints
            con = [obj.A - eye(size_x) obj.B;obj.C zeros(size_y,size_u)]*[S_c;R_c]==[obj.P*obj.B_hat;zeros(size_y,size_uhat)];
            
            % solve matrices S and R considering equation (17)
            ops = sdpsettings('verbose',obj.opt_option.verbose);
            sol = optimize(con,'', ops);
            
            if sol.problem == 0
                % solve sucessfully considering equation (17)
                fea = 1;
                obj.S = value(S_c);
                obj.R = value(R_c);
            else
                % solve equation (17) uncessfully, use (11) instead
                fea = 0;
                B_bar = [obj.A-eye(size_x) obj.B];
                SR = (B_bar'*obj.M*B_bar)\(B_bar'*obj.M*obj.P*obj.B_hat);% solving (11) similar to (16)  
                obj.S = SR(1:size_x,:);
                obj.R = SR(size_x+1:size_x+size_u,:);
            end
        end
        
        %% checking the feasibility of the relation leveraging equation (15)
        function [obj,fea] = check_fea(obj)
            % check_fea: compute gamma1 in (11) and gamma2 in (12) and check the feasibility condition in (15)
            % Input:    obj: a gAAS_relation_options object
            % output:   fea: flag indicating whether the feasibility condition in (15) holds 
            %               fea = 1: feasibility condition is respected
            %               fea = 0: feasibility condition is NOT respected
            %               fea = -1: problem occurs when computing gamma1 or gamma2 
            %           obj: a gAAS_relation_options object with the computed gamma1 and gamma2
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % compute gamma1
            gamma1_cur = obj.compute_gamma1();
            
            % compute gamma2
            gamma2_cur = obj.compute_gamma2();
            
            if isempty(gamma1_cur) || isempty(gamma2_cur)
                % gamma1 or/gamma2 cannot be computed due to unknown reason
                fea = -1;
                disp('Checking feasibility condition: gamma1 or/and gamma2 cannot be computed, please check the problem setting');
            else
                obj.gamma1 = gamma1_cur;
                obj.gamma2 = gamma2_cur;
                if gamma1_cur + gamma2_cur <= obj.varepsilon*(1-sqrt(obj.kappa))
                    % feasibility condition in equation (15)
                    fea = 1;
                    disp('Checking feasibility condition: feasibility condition for gamma1 and gamma2 is respected.');
                else
                    fea = 0;
                    disp('Checking feasibility condition: feasibility condition for gamma1 and gamma2 is NOT respected.');
                end
                
            end
        end
        
        % compute gamma1 in equation (11)
        function gamma1 = compute_gamma1(obj)
            % compute_gamma1: compute gamma1 in equation (11)
            % Input:    obj: a gAAS_relation_options object
            % output:   gamma1: gamma1 as in equation (11)
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % define variable for uhat
            size_uhat = size(obj.B_hat,2);
            u_hat = sdpvar(size_uhat,1);
            
            % adding all input constraints to the optimization problem
            num_Uhat = size(obj.U_hat);
            con = [];
            for i = 1:1:num_Uhat
                con = [con;obj.U_hat{1}.A*u_hat<=obj.U_hat{1}.b];
            end
            
            % specifying the obtimization problem for solving gamma1
            objfun = -u_hat'*(obj.A*obj.S+obj.B*obj.R-obj.P*obj.B_hat-obj.S)'*obj.M*(obj.A*obj.S+obj.B*obj.R-obj.P*obj.B_hat-obj.S)*u_hat;
            ops = sdpsettings('solver', 'fmincon', 'verbose',obj.opt_option.verbose);
            sol = optimize(con, objfun, ops);
            
            if sol.problem == 0
                % solve sucessfully
                gamma1 = -value(objfun);
            else
                % solve unsucessfully due to unknow reason
                gamma1 = [];
            end
        end
        
        % compute gamma2 in equation (12)
        function gamma2 = compute_gamma2(obj)
            % compute_gamma2: compute gamma2 in equation (12)
            % Input:    obj: a gAAS_relation_options object
            % output:   gamma2: gamma1 as in equation (12)
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % define variables for uhat
            size_xhat = size(obj.A_hat,2);
            x_hat = sdpvar(size_xhat,1);
            
            % adding all state constraints to the optimization problem
            num_Xhat = size(obj.X_hat);
            con = [];
            for i = 1:1:num_Xhat
                con = [con;obj.X_hat{1}.A*x_hat<=obj.X_hat{1}.b];
            end
            
            % specifying the obtimization problem for solving gamma2
            objfun = -x_hat'*obj.D'*obj.M*obj.D*x_hat;
            ops = sdpsettings('solver', 'fmincon', 'verbose',obj.opt_option.verbose);
            sol = optimize(con, objfun, ops);
            
            if sol.problem == 0
                % solve sucessfully
                gamma2 = -value(objfun);
            else
                % solve unsucessfully due to unknow reason
                gamma2 = [];
            end
        end
        
        %% Computing control input for the abstract system considering Lemma 4.1 and Remark 4.3
        function u_hat_cur = provide_u(obj,u_prevs,u_sug)
            % provide_u: providing control inputs for the abstract system that respecting the constraints in Lemma 4.1 to maintain the gAAS relation
            % Input:    obj: a gAAS_relation_options object
            % output:   u_hat_cur: control inputs for the abstract system that respect the gAAS relation
            %   Authors: Bingzhuo Zhong
            %   Date: April 1, 2024
            
            % first give a fast check whether the current suggested control
            % input falls in the constraint as described in Lemma 4.1
            if (u_sug-u_prevs)'*obj.S'*obj.M*obj.S*(u_sug-u_prevs)<= (obj.varepsilon*(1-sqrt(obj.kappa))-obj.gamma1-obj.gamma2)^2
                u_hat_cur = u_sug;
            else 
                % if not, solve a control input that is the closest to the suggested one and respecting the constraints in Lemma 4.1
                u = sdpvar(size(obj.B,2),1);
                optfun = (u-u_sug)'*(u-u_sug);
                con= (u-u_prevs)'*obj.S'*obj.M*obj.S*(u-u_prevs)<= (obj.varepsilon*(1-sqrt(obj.kappa))-obj.gamma1-obj.gamma2)^2;
                ops = sdpsettings('solver', 'fmincon', 'verbose',obj.opt_option.verbose);
                sol = optimize(con, optfun, ops);
                                
            if sol.problem == 0
                % solve sucessfully
                u_hat_cur = value(u);
            else
                % if error occurs in the solution process (very unlikely), use the previous input which definitely satisfies the constraints specified in Lemma 4.1
                u_hat_cur = u_prevs;
            end
            end            
        end
        
    end
end

