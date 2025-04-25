function [M_ParamMarginal_SimulationPoints,M_ParamCovariance_SimulationPoints]=Interp_ModelParameters(V_X_ObsSites,V_Y_ObsSites,M_ParamMarginal_ObsSites,M_ParamCovariance_ObsSites,V_X_SimulationPoints,V_Y_SimulationPoints)
%******
%Interpolates the parameters of the stochastic rainfall model at simulation locations from observation sites
%******

nb_clusters=size(M_ParamMarginal_ObsSites,3);

%Step 3.1: Interpolate marginal model on SimulationGrid (Ordinary Kriging)
M_ParamMarginal_SimulationPoints=NaN(length(V_X_SimulationPoints),nb_clusters,3); %also simulate at observation locations (for later conditioning)
str_paramTPS_marginal=struct();
for my_type=1:nb_clusters
    for my_param=1:3
        %set constraints
        if my_param==1
            str_paramTPS_marginal(my_type,my_param).min_constr=-4;
            str_paramTPS_marginal(my_type,my_param).max_constr=4;
        elseif my_param==2
            str_paramTPS_marginal(my_type,my_param).min_constr=0.01;
            str_paramTPS_marginal(my_type,my_param).max_constr=25;
        else
            str_paramTPS_marginal(my_type,my_param).min_constr=1.01;
            str_paramTPS_marginal(my_type,my_param).max_constr=250;
        end
        %Ordinary Krigin (Matérn covariance with nugget)
        [V_param]=OK_interp_MarginalParams(V_X_ObsSites, V_Y_ObsSites, M_ParamMarginal_ObsSites(:,my_param,my_type), V_X_SimulationPoints, V_Y_SimulationPoints);
        str_paramTPS_marginal(my_type,my_param).alpha_TPS=-1;
        M_ParamMarginal_SimulationPoints(:,my_type,my_param)=V_param;
    end
end

%Step 3.2: Interpolate covariance model on SimulationGrid (Thin Plate Interpolation)

M_ParamCovariance_SimulationPoints=zeros(length(V_X_SimulationPoints),nb_clusters,5);
str_paramTPS_covariance=struct();
for my_type=1:nb_clusters
    for my_param=1:5

        %TPS interpolation
        if my_param<5
            if my_param==1
                str_paramTPS_covariance(my_type,my_param).min_constr=-150;
                str_paramTPS_covariance(my_type,my_param).max_constr=150;
            elseif my_param==2
                str_paramTPS_covariance(my_type,my_param).min_constr=-150;
                str_paramTPS_covariance(my_type,my_param).max_constr=150;
            elseif my_param==3
                str_paramTPS_covariance(my_type,my_param).min_constr=5;
                str_paramTPS_covariance(my_type,my_param).max_constr=150;
            else
                str_paramTPS_covariance(my_type,my_param).min_constr=0.3;
                str_paramTPS_covariance(my_type,my_param).max_constr=2;
            end
            [V_param,alpha_TPS]=TPS_Interp(V_X_ObsSites,V_Y_ObsSites,M_ParamCovariance_ObsSites(:,my_type,my_param),V_X_SimulationPoints,V_Y_SimulationPoints,str_paramTPS_covariance(my_type,my_param).min_constr,str_paramTPS_covariance(my_type,my_param).max_constr,0.01,0);
            str_paramTPS_covariance(my_type,my_param).alpha_TPS=alpha_TPS;
            M_ParamCovariance_SimulationPoints(:,my_type,my_param)=V_param;
        end

    end
end

end