function[M_iter_simul_latent]=StochasticSimulation_NonStationaryGaussianField(my_type,V_X_SimulationGrid,V_Y_SimulationGrid,M_ParamCovariance_SimulationGrid,nb_realiz)

%Init
M_iter_simul_latent=NaN(length(V_X_SimulationGrid),nb_realiz); %lines: simulation locations, columns: realizations

%Choleski factorization of the covariance matrix on simulationGrid
V_g1_target=M_ParamCovariance_SimulationGrid(:,my_type,1);
V_g2_target=M_ParamCovariance_SimulationGrid(:,my_type,2);
V_rho2_target=M_ParamCovariance_SimulationGrid(:,my_type,3);
V_nu_target=M_ParamCovariance_SimulationGrid(:,my_type,4);

Sigma_TT=NonStationaryAnisotropic_MaternCov(V_X_SimulationGrid/111,V_X_SimulationGrid/111,V_Y_SimulationGrid/111,V_Y_SimulationGrid/111,V_g1_target,V_g1_target,V_g2_target,V_g2_target,V_rho2_target,V_rho2_target,V_nu_target,V_nu_target);
disp('Choleski factorization')
det(Sigma_TT)
L=chol(Sigma_TT,'lower');

disp('Simulation')
for my_sim=1:nb_realiz
    Z_Sim=L*randn(length(V_X_SimulationGrid),1);
    M_iter_simul_latent(:,my_sim)=Z_Sim;
end

end