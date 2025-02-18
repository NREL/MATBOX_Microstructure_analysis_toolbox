%% FRACTAL DIMENSION (COUNTING BOX)
if p.fractal_boxcounting.todo
    fprintf('> Phase domain fractal dimension (box-counting method)\n');
    p.fractal_boxcounting.topology_dimension = number_dimension;
    p.fractal_boxcounting.plot = false;
    for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
        binary_array = zeros(size(Phase_microstructure));
        binary_array(Phase_microstructure==phaselabel(current_phase_todo,1))=1;
        [N(:,current_phase_todo),box_lengths,fractal_dimension(:,current_phase_todo),fractal_dimension_convergence(:,:,current_phase_todo)] = Function_fractaldimension_boxcounting(binary_array,p.fractal_boxcounting);
    end

    % Table
    Table_fractaldimension = table(phasename_todo,fractal_dimension(1,:)',fractal_dimension(2,:)',fractal_dimension(3,:)',fractal_dimension(4,:)',...
        'VariableNames',{'Phase' 'Fit from 1 to' 'Fractal dimension' 'Topology dimension' 'Fractal propensity'});
    for current_phase_todo=1:1:number_phase_todo
        Table_boxlength(current_phase_todo).t = table(fractal_dimension_convergence(:,1,current_phase_todo),fractal_dimension_convergence(:,2,current_phase_todo),fractal_dimension_convergence(:,3,current_phase_todo),...
        'VariableNames',{'Fit from 1 to' 'Fractal dimension' 'Fit norm error'});
    end
    Results_Volumefraction.Table_fractaldimension = Table_fractaldimension; % Save in main table result
    disp(Table_fractaldimension)

    % Save table
    if opt.save.xls
        filename = 'PhaseDomain_Fractaldimension_boxcounting'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Fractal dimension';
        DATA_writetable.sheet(1).table = Table_fractaldimension;
        for current_phase_todo=1:1:number_phase_todo
            DATA_writetable.sheet(1+current_phase_todo).name = char(phasename_todo(current_phase_todo));
            DATA_writetable.sheet(1+current_phase_todo).table = Table_boxlength(current_phase_todo).t;
        end           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    % Figure
    if p.fractal_boxcounting.topology_dimension==3
        propertyname= 'Phase volume'; % Figure name
    else
        propertyname = 'Phase surface'; % Figure name
    end
    filename = 'PhaseDomain_fractal_dimension_boxcounting';
    function_fractalfig(propertyname,filename,Current_folder, box_lengths,N,fractal_dimension_convergence,p.fractal_boxcounting.topology_dimension,number_phase,p,opt,infovol);

    % Correlation
    results_correlation(current_phase_todo).PhaseDomain_fractaldimension_boxcounting = fractal_dimension(2,current_phase_todo);
    results_correlation(current_phase_todo).PhaseDomain_fractalpropensity_boxcounting = abs(p.fractal_boxcounting.topology_dimension - fractal_dimension(2,current_phase_todo));
end