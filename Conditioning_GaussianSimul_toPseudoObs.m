function[M_simul_condi_latent]=Conditioning_GaussianSimul_toPseudoObs(my_type,M_ToCondition_latent,V_X_SimilationGrid,V_Y_SimilationGrid,M_ParamCovariance_SimulationGrid,V_X_Condi,V_Y_Condi,M_ParamCovariance_Condi,V_condi_data_latent)

disp('start preliminary processing')
%Preliminary1: Covariance matrix of observations (and its inverse)
[Sigma]=NonStationaryAnisotropic_MaternCov(V_X_Condi/111,V_X_Condi/111,V_Y_Condi/111,V_Y_Condi/111,M_ParamCovariance_Condi(:,my_type,1),M_ParamCovariance_Condi(:,my_type,1),M_ParamCovariance_Condi(:,my_type,2),M_ParamCovariance_Condi(:,my_type,2),M_ParamCovariance_Condi(:,my_type,3),M_ParamCovariance_Condi(:,my_type,3),M_ParamCovariance_Condi(1,my_type,4),M_ParamCovariance_Condi(1,my_type,4));
Cab=Sigma;
inv_cab=inv(Cab);

%Preliminary2: Covariance matrix between observations and simulations
[Cav]=NonStationaryAnisotropic_MaternCov(V_X_Condi/111,V_X_SimilationGrid/111,V_Y_Condi/111,V_Y_SimilationGrid/111,M_ParamCovariance_Condi(:,my_type,1),M_ParamCovariance_SimulationGrid(:,my_type,1),M_ParamCovariance_Condi(:,my_type,2),M_ParamCovariance_SimulationGrid(:,my_type,2),M_ParamCovariance_Condi(:,my_type,3),M_ParamCovariance_SimulationGrid(:,my_type,3),M_ParamCovariance_Condi(1,my_type,4),M_ParamCovariance_SimulationGrid(1,my_type,4));


%Preliminary3: Assign corresponding grid cell to each observation site
V_ind_correspondingGridCell=(length(V_X_SimilationGrid)-length(V_X_Condi)+1):length(V_X_SimilationGrid); %conditioning location are included in unconditional simulations at the end of the SimulationGrid

disp('end preliminary processing')

M_simul_condi_latent=NaN(size(M_ToCondition_latent));
%Conditioning Kriging
for my_realiz=1:size(M_ToCondition_latent,2)
    %Simple Kriging obs
    Za=V_condi_data_latent(:,min(my_realiz,size(V_condi_data_latent,2)));
    Lambda_a=inv_cab*Cav;
    V_Kriged_Obs=Za'*Lambda_a;
    %Simple Kriging uncondi simulation
    V_uncondi_SiteObs=M_ToCondition_latent(V_ind_correspondingGridCell,my_realiz);
    Za=V_uncondi_SiteObs;
    Lambda_a=inv_cab*Cav;
    V_Kriged_Sim=Za'*Lambda_a;
    %Perform conditioning
    M_simul_condi_latent(:,my_realiz)=V_Kriged_Obs'+M_ToCondition_latent(:,my_realiz)-V_Kriged_Sim';
end
disp('end conditioning')

end