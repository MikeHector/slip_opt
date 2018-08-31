%Let's try to make this analytic
clear; close all; clc;

makefile = 1; makehessian = 0;
%Set up collocation parameters - these are fixed for every compile
p.N = 90;
p.Nstance = 60;
p.Nflight = 30;
p.dof = 3;
p.cntrl_dof = 2;

% Setup
[dv, vStruc, vList, Example] = makeParamSym(p);

% Objective function
%Get symbolic: objective func - gradient - hessian
obj_func = objectiveFunctionSymbolic(dv, p, vStruc);
grad_obj = jacobian(obj_func, dv).';
hess_obj = cell(1,size(grad_obj,2));
for i = 1:size(grad_obj,2)
    hess_obj{i} = jacobian(grad_obj(:,i),dv);
end

if makefile == 1
    %Objective function file
    currdir = [pwd filesep];
    filename = [currdir, 'SymObjFunc.m'];
    disp('Making objective + grad file')
    matlabFunction(obj_func,grad_obj,'file',filename,'vars',{dv, vList});
    disp('Done')
    
    if makehessian == 1
        %Objective Hessian
        filename = [currdir, 'SymObjHessFunc.m'];
        disp('making objective hessian file')
        matlabFunction(hess_obj{size(grad_obj,2)},'file',filename,'vars',{dv, vList});
        disp('Done')
    end
end

% Constraints
%Get symbolic constraints
[c_ineq, c_eq, c_linear, A, b] = getSymbolicConstraints(dv,p,vStruc);

% [c_ineq_from_old, c_eq_from_old] = nonlinear_constraint_func2(dv_num_block, symparamEX);

if makefile == 1
    filename = 'c_ineq_func';
    disp('Making c_ineq file')
    matlabFunction(c_ineq,'file',filename,'vars',{dv, vList},...
        'outputs',{'c_ineq'}, 'Optimize', false);
    disp('Done')
end

if makefile == 1
    filename = 'c_eq_func';
    disp('Making c_eq file')
    matlabFunction(c_eq,'file',filename,'vars',{dv, vList},...
        'outputs',{'c_eq'},'Optimize', false);
end

%Now we compare!

%Get old 
[c_ineq_old_block, c_eq_old_block] = nonlinear_constraint_func2(Example.dvNum, Example.vNum);
%Get it into list form
c_ineq_old_list = [c_ineq_old_block(1,:)'; c_ineq_old_block(2,:)'];
c_eq_old_list = []; %zeros( size(c_eq_old_block,1)*size(c_eq_old_block,2) ,1);
for i = 1:size(c_eq_old_block,2)
    c_eq_old_list = [c_eq_old_list; c_eq_old_block(:,i)];
end

%Get new
c_ineq_new = c_ineq_func(Example.dvSym, Example.vSymList);
c_eq_new = c_eq_func(Example.dvSym, Example.vSymList);



disp('stop'); pause; pause; pause(60); pause;
grad_c_eq = jacobian(c_eq, DV).';
grad_c_ineq = jacobian(c_ineq,DV).';

%Constraint function files
filename = [currdir, 'SymConFunc.m'];
matlabFunction(c_ineq,c_eq,grad_c_ineq,grad_c_eq,'file',filename,'vars',{DV, symparamlist},...
    'outputs',{'c_ineq','c_eq','grad_c_ineq','grad_c_eq'});

i_c_eq = size(grad_c_eq,2);
hess_c_eq = cell(1,i_c_eq);
for i = 1:i_c_eq
    hess_c_eq{i} = jacobian(grad_c_eq(:,i),DV);
    disp('Constraint Hessian 1 of 2 generating...')
    disp(['Progress: ', num2str(i/i_c_eq*100), ' %'])
end

i_c_ineq = size(grad_c_ineq,2);
hess_c_ineq = cell(1,i_c_ineq);
for i = 1:i_c_ineq
    hess_c_ineq{i} = jacobian(grad_c_ineq(:,i),DV);
    disp('Constraint Hessian 2 of 2 generating...')
    disp(['Progress: ', num2str(i/i_c_ineq*100), ' %'])
end

%%%
%Constraint Hessian

% %% Build files

% %%%

time2gen = toc;
disp(['It took ' num2str(time2gen/60), ' minutes to generate analytics'])


