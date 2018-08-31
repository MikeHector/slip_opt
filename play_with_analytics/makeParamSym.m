function [dv, PsymStruc, v, Example] = makeParamSym(CollocationParameters)
%makeParamSym Make inputs for creating analytical derivatives, comparing
%them

    %Make the decision variable symbols as a long list
    numberOfCollocationPoints = (2* CollocationParameters.dof +...
        CollocationParameters.cntrl_dof) * (CollocationParameters.N -1) + 2;
    dv = sym('x', [1 numberOfCollocationPoints], 'real')'; %all my decision vars (one hop with flight)
    %x1 = x(1), x2 = y(1) ... x8 = Tankle(1) x9 = x(2) ... ... ... ... ... ...
    %x(end -1) = Tstance x(end) = Tflight

    %Make the parameters that will be passed as symbols as well
    syms R_leg R_ankle c m k i_motor transmission lf transmission_ankle g disturbance_f apex_height apex_velocity real
    v = [m c k g i_motor transmission transmission_ankle lf ...
                R_leg R_ankle disturbance_f apex_height apex_velocity]';
    v_list = [30 10 2000 10 .00005 80 50 .15 1 3 3 1 1]'; %Example parameter list

    %Let Psymparam inheirit properties of CollocationParameters
    PsymStruc = CollocationParameters;
    PsymStrucEX = CollocationParameters;
    %Compile symbolic and example array into struct
    for i=1:length(v)
        PsymStruc.(char(v(i))) = v(i);
        PsymStrucEX.(char(v(i))) = v_list(i);
    end
    
    %Random set of decision vars for checking later
    dv_num_list = rand([length(dv),1]);
    rowLen = 2*CollocationParameters.dof...
            +CollocationParameters.cntrl_dof;
    %Build into block
    for i = 1:CollocationParameters.N-1
        dv_num_block(:,i) = dv_num_list((i-1) * rowLen+1: (i-1) * rowLen+8);
    end
    dv_num_block(end+1,1:2) = [dv_num_list(end-1), dv_num_list(end)];
    
    Example.dvSym = dv_num_list;
    Example.vSymList = v_list;
    Example.dvNum = dv_num_block;
    Example.vNum = PsymStrucEX;
    
end