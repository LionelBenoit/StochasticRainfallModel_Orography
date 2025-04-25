function [V_X_MonthlyCondi,V_Y_MonthlyCondi]=Set_ConditioningLocations_Monthly(V_X_Condi,V_Y_Condi,V_X_SimulationGrid,V_Y_SimulationGrid,nb_condiMonthly)
    
V_X_Condi_tot=V_X_Condi;
V_Y_Condi_tot=V_Y_Condi;

for i=1:nb_condiMonthly
    M_dist=pdist2([V_X_SimulationGrid V_Y_SimulationGrid],[V_X_Condi_tot V_Y_Condi_tot]);
    V_max=min(M_dist,[],2);
    my_ind=find(V_max==max(V_max));
    my_ind=my_ind(1);
    V_X_Condi_tot=[V_X_Condi_tot;V_X_SimulationGrid(my_ind)];
    V_Y_Condi_tot=[V_Y_Condi_tot;V_Y_SimulationGrid(my_ind)];
end

V_X_MonthlyCondi=V_X_Condi_tot(length(V_X_Condi)+1:end);
V_Y_MonthlyCondi=V_Y_Condi_tot(length(V_Y_Condi)+1:end);

end