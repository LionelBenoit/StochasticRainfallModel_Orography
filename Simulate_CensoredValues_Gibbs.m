function[M_Data_latent_full]=Simulate_CensoredValues_Gibbs(M_Data_latent_censored,V_X,V_Y,M_ParamCovariance,M_ParamMarginal,V_RainTypes)

M_Data_latent_full=M_Data_latent_censored;

for my_type=1:max(V_RainTypes)
    my_type
    [Sigma]=NonStationaryAnisotropic_MaternCov(V_X/111,V_X/111,V_Y/111,V_Y/111,M_ParamCovariance(:,my_type,1),M_ParamCovariance(:,my_type,1),M_ParamCovariance(:,my_type,2),M_ParamCovariance(:,my_type,2),M_ParamCovariance(:,my_type,3),M_ParamCovariance(:,my_type,3),M_ParamCovariance(:,my_type,4),M_ParamCovariance(:,my_type,4));
    inv_Sigma_Gibbs=inv(Sigma);
    Cab=Sigma;
    inv_cab=inv(Cab);

    inds_type=find(V_RainTypes==my_type);

    %inequality constraint Gibbs sampling to simulate censored values
    for my_day=1:length(inds_type)
        display(strcat(num2str(my_day),'/',num2str(length(inds_type))))
        %my_day
        %initialization
        V_data=M_Data_latent_censored(:,inds_type(my_day));
        V_data_simulGibbs=V_data;
        for my_pt=1:length(V_data)
            if V_data_simulGibbs(my_pt)==-5
                V_data_simulGibbs(my_pt)=M_ParamMarginal(my_pt,1,my_type)-0.1;
            end
        end
        %Gibbs sampling
        for my_iter=1:2500
            for my_pt=1:length(V_data)
                if V_data(my_pt)==-5
                    [ mu_gibbs, sigma_gibbs ] = conditional_normal_sim(inv_Sigma_Gibbs,my_pt,V_data_simulGibbs);
                    sim_value=mu_gibbs+randn*sigma_gibbs;
                    it=1;
                    while sim_value > M_ParamMarginal(my_pt,1,my_type)

                        sim_value=mu_gibbs+randn*sigma_gibbs;
                        it=it+1;
                        if it>100
                            sim_value=V_data_simulGibbs(my_pt);
                            break
                        end
                    end
                    V_data_simulGibbs(my_pt)=sim_value;
                end
            end
        end
        M_Data_latent_full(:,inds_type(my_day))=V_data_simulGibbs;
    end

end

end %end function