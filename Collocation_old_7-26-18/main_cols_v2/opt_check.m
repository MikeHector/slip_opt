%Optimality check
%For each optimal traj check a couple different seeds to make sure they
%converge to the same solution

strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');

for i = 1:length(strucc)
    filename = strucc(i).name;
    filename = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\', filename);
    load(filename)
    results{i} = opt_results;
    c(i) = opt_results.c;
end
[c_sorted,i] = sort(c);

q=1;

for k = 1:length(i)
    results_sorted_c{k} = results{i(k)};
    flags(k) = results{i(k)}.flag;
    if results{i(k)}.flag > 0
        c_graph(q) = results{i(k)}.c;
        cost_graph(q) = results{i(k)}.cost;
        q = q+1;
    end
end

[~, keep] = spaced(flags, c_sorted);
BADTHINGSARECOMING = 0;
for i = 257:numel(results_sorted_c)
    if ismember(results_sorted_c{i}.c, keep)
        %Check some seeds
        legacyCost = results_sorted_c{i}.cost;
        X = results_sorted_c{i}.X;
        res = results_sorted_c{i};
        seed1 = zeros(size(X));
        seed2 = X + rand(size(X));
        seed3 = rand(size(X));
        seed4 = [linspace(-.2, .2, length(X));...%x
                 0 * linspace(1, 2, length(X));...%y
                 linspace(.9, 1.1, length(X));...%r0
                 linspace(1, 1, length(X));...%dx
                 0 * linspace(1, 2, length(X));...%dy
                 0 * linspace(1, 2, length(X));...%dr0
                 0 * linspace(1, 2, length(X));...%Tleg
                 0 * linspace(1, 2, length(X));...%Tankle
                 [.2, linspace(0,0, length(X)-1)]];%Time
         seeds = {seed2, seed3, seed4}; %{seed1, seed2, seed3, seed4};
         for k = 1:numel(seeds)
            X0 = seeds{k};
            opt_results.flag = 0;
            counter = 0;
            while opt_results.flag == 0 && counter < 10
                [X0, opt_results] = RUN_COL(X0, res.c, res.apex_velocity, res.apex_height, 1, res.apex_velocity, 0, atan2(X(2,1),X(1,1)) - pi/2);
                counter = counter + 1;
            end
            filename = strcat('optCheck_damping_', num2str(cputime*10000000));
            save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\',filename),'opt_results');
            if opt_results.cost < legacyCost     
                BADTHINGSARECOMING = BADTHINGSARECOMING + 1;
            end
         end
    end
end