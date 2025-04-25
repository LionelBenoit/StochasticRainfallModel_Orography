function[M_Data_latent_MwG_Results,M_RainValues_MwG_Results]=Metropolis_within_Gibbs(inds_target_month,M_DataCondi_Month,V_RainTypes_Month,V_X_PseudoObsSites,V_Y_PseudoObsSites,V_Obs_MonthlyCondi,M_ParamMarginal_PseudoObs,M_ParamCovariance_PseudoObs,V_X_ObsSites,V_Y_ObsSites,M_ParamCovariance_ObsSites,nb_realiz_MwG,nb_iter_MwG)

nb_days=length(inds_target_month);

%Extract data of the target month
V_RandomPath_days=randperm(length(inds_target_month));

%Pre-compute inv_Sigma to save computation time
nb_clusters=max(V_RainTypes_Month);
M_invSigma_MwG=NaN(length(V_X_PseudoObsSites),length(V_X_PseudoObsSites),nb_clusters);
for my_type_init=1:nb_clusters
    [Sigma]=NonStationaryAnisotropic_MaternCov(V_X_PseudoObsSites/111,V_X_PseudoObsSites/111,V_Y_PseudoObsSites/111,V_Y_PseudoObsSites/111,M_ParamCovariance_PseudoObs(:,my_type_init,1),M_ParamCovariance_PseudoObs(:,my_type_init,1),M_ParamCovariance_PseudoObs(:,my_type_init,2),M_ParamCovariance_PseudoObs(:,my_type_init,2),M_ParamCovariance_PseudoObs(:,my_type_init,3),M_ParamCovariance_PseudoObs(:,my_type_init,3),M_ParamCovariance_PseudoObs(:,my_type_init,4),M_ParamCovariance_PseudoObs(:,my_type_init,4));
    inv_Sigma=inv(Sigma);
    M_invSigma_MwG(:,:,my_type_init)=inv_Sigma;
end

%Pre-compute unconditional simulations
M_SimulUncondi_PseudoObs_AllTypes=NaN(length(V_X_PseudoObsSites),nb_realiz_MwG*nb_days,nb_clusters);
for my_type=1:nb_clusters
    [M_SimulUncondi_PseudoObs]=StochasticSimulation_NonStationaryGaussianField(my_type,V_X_PseudoObsSites,V_Y_PseudoObsSites,M_ParamCovariance_PseudoObs,nb_realiz_MwG*nb_days);
    M_SimulUncondi_PseudoObs_AllTypes(:,:,my_type)=M_SimulUncondi_PseudoObs;
end

%Initialize MwG matrices

M_Data_latent_MwG=[zeros(length(V_X_PseudoObsSites)-length(V_X_ObsSites),nb_days); M_DataCondi_Month];
M_RainValues_MwG=zeros(length(V_X_PseudoObsSites),nb_days);

M_Data_latent_MwG_Results=NaN([size(M_Data_latent_MwG) nb_realiz_MwG]);
M_RainValues_MwG_Results=NaN([size(M_RainValues_MwG) nb_realiz_MwG]);

for my_realiz_MwG=1:nb_realiz_MwG
    my_realiz_MwG

    %Initialization MwG by conditional simulation
    for my_day=1:nb_days
        my_day
        %unconditional simulation
        my_type=V_RainTypes_Month(my_day);
        if my_type>0

            %unconditional simulation
            V_SimulUncondi_PseudoObs_currentDay=M_SimulUncondi_PseudoObs_AllTypes(:,my_day+(my_realiz_MwG-1)*nb_days,my_type);
            %[V_SimulUncondi_PseudoObs_currentDay]=StochasticSimulation_NonStationaryGaussianField(my_type,V_X_PseudoObsSites,V_Y_PseudoObsSites,M_ParamCovariance_PseudoObs,1);

            %conditioning to observations (i.e. in the latent world)
            V_condi_data_latent=M_DataCondi_Month(:,my_day); %!!censored values in M_DataCondi_latent must have been preliminarily simulated by Gibbs sampling
            [V_simul_condi_latent]=Conditioning_GaussianSimul_toPseudoObs(my_type,V_SimulUncondi_PseudoObs_currentDay,V_X_PseudoObsSites,V_Y_PseudoObsSites,M_ParamCovariance_PseudoObs,V_X_ObsSites,V_Y_ObsSites,M_ParamCovariance_ObsSites,V_condi_data_latent);
            M_Data_latent_MwG(:,my_day)=V_simul_condi_latent;

            %Back-transform
            V_Y=V_simul_condi_latent;
            [V_R]=Transform_Y_to_R(V_Y,M_ParamMarginal_PseudoObs(:,my_type,1),M_ParamMarginal_PseudoObs(:,my_type,2),M_ParamMarginal_PseudoObs(:,my_type,3));
            M_RainValues_MwG(:,my_day)=V_R;
        end
    end

    for my_pt=1:(length(V_X_PseudoObsSites)-length(V_X_ObsSites))
        val_monthlyMap=V_Obs_MonthlyCondi(my_pt);

        V_Rain_current=M_RainValues_MwG(my_pt,:);
        val_dailyMaps=sum(V_Rain_current);

        if val_monthlyMap>1 && val_dailyMaps>1
            coef=val_monthlyMap/val_dailyMaps;
            V_Rain_current=V_Rain_current*coef;
            M_RainValues_MwG(my_pt,:)=V_Rain_current;
            for my_day=1:nb_days
                my_type=V_RainTypes_Month(my_day);
                if (V_Rain_current(my_day)>0) && (my_type>0)
                    my_latent_val=Transform_R_to_Y(V_Rain_current(my_day),M_ParamMarginal_PseudoObs(my_pt,my_type,1),M_ParamMarginal_PseudoObs(my_pt,my_type,2),M_ParamMarginal_PseudoObs(my_pt,my_type,3));
                    M_Data_latent_MwG(my_pt,my_day)=my_latent_val;
                end
            end
        end
    end

    %MwG per se
    nb_tests=0;
    nb_accept=0;
    for my_iter=1:nb_iter_MwG
        if mod(my_iter,50)==0
            display(strcat('Realiz=',num2str(my_realiz_MwG),'-Iter=',num2str(my_iter),'-Accept=',num2str(nb_accept/nb_tests*100),'%'))
        end
        %Gibbs sampling
        for my_day_loop=1:nb_days

            my_day=V_RandomPath_days(my_day_loop); %explore days following a random path
            my_type=V_RainTypes_Month(my_day);

            if my_type>0

                inv_Sigma_Gibbs=M_invSigma_MwG(:,:,my_type);

                for my_pt=1:(length(V_X_PseudoObsSites)-length(V_X_ObsSites))
                    [ mu_gibbs, sigma_gibbs ] = conditional_normal_sim(inv_Sigma_Gibbs,my_pt,M_Data_latent_MwG(:,my_day));
                    sim_value_Gibbs=mu_gibbs+randn*sigma_gibbs;

                    V_Rain_current=M_RainValues_MwG(my_pt,:);

                    V_Rain_test=V_Rain_current;
                    rainValue_test=Transform_Y_to_R(sim_value_Gibbs,M_ParamMarginal_PseudoObs(my_pt,my_type,1),M_ParamMarginal_PseudoObs(my_pt,my_type,2),M_ParamMarginal_PseudoObs(my_pt,my_type,3));
                    V_Rain_test(my_day)=rainValue_test;


                    Monthly_sum_current=sum(V_Rain_current);
                    Monthly_sum_test=sum(V_Rain_test);

                    %Metropolis rule
                    my_sig=min(2,V_Obs_MonthlyCondi(my_pt)*0.005);
                    %my_sig=10;

                    L0=1/(my_sig*sqrt(2*pi))*exp(-1/(2*my_sig*my_sig)*(Monthly_sum_current-V_Obs_MonthlyCondi(my_pt))^2);
                    L1=1/(my_sig*sqrt(2*pi))*exp(-1/(2*my_sig*my_sig)*(Monthly_sum_test-V_Obs_MonthlyCondi(my_pt))^2);

                    alpha=L1/L0;
                    nb_tests=nb_tests+1;

                    if alpha>1 || rand > 1-alpha
                        %display('accept')
                        M_Data_latent_MwG(my_pt,my_day)=sim_value_Gibbs;
                        M_RainValues_MwG(my_pt,my_day)=rainValue_test;
                        nb_accept=nb_accept+1;
                    end

                end

            end
        end

    end
    M_Data_latent_MwG_Results(:,:,my_realiz_MwG)=M_Data_latent_MwG;
    M_RainValues_MwG_Results(:,:,my_realiz_MwG)=M_RainValues_MwG;

end

end