%Sub-script to be called in Script_Dailymaps_NonStationaryMarginalAndCov.m

for my_type=1:nb_clusters
%for my_type=2:2

    M_Data_type=M_Data_Training(:,V_RainTypes_Training==my_type);

    figure(num_fig+my_type)
    clf
    colormap('jet')
    
    
    %Marginal - a0 - at simulation locations
    subplot(2,7,1)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_gauge=1:size(M_Data_gapfilled,1)
        my_max=6;
        my_val=M_ParamMarginal_ObsSites(my_gauge,1,my_type)+3;
        my_ind_color=max(1,min(floor(my_val/my_max*255)+1,255));
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([-3 3])
    title('a0')

    %Marginal - k - at simulation locations
    subplot(2,7,2)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_gauge=1:size(M_Data_gapfilled,1)
        my_max=1.5;
        my_val=M_ParamMarginal_ObsSites(my_gauge,2,my_type);
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('k')

    %Marginal - teta - at simulation locations
    subplot(2,7,3)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_gauge=1:size(M_Data_gapfilled,1)
        my_max=25;
        my_val=M_ParamMarginal_ObsSites(my_gauge,3,my_type);
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('teta')

    %Covariance - gamma1 - at PilotPoint locations
    subplot(2,7,4)
    hold on
    my_max=150;
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_pilotpoint=1:length(V_X_ClimateDivisions)
        my_val=M_ParamsMatern_ClimateDivisions(my_type,my_pilotpoint,1)+150;
        my_ind_color=max(min(floor(my_val/300*255)+1,255),1);
        plotm(V_Y_ClimateDivisions(my_pilotpoint)/111,V_X_ClimateDivisions(my_pilotpoint)/111,'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([-my_max my_max])
    title('gamma1')
    

    %Covariance - gamma2 - at PilotPoint locations
    subplot(2,7,5)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_pilotpoint=1:length(V_X_ClimateDivisions)
        my_val=M_ParamsMatern_ClimateDivisions(my_type,my_pilotpoint,2)+150;
        my_ind_color=max(min(floor(my_val/300*255)+1,255),1);
        plotm(V_Y_ClimateDivisions(my_pilotpoint)/111,V_X_ClimateDivisions(my_pilotpoint)/111,'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([-my_max my_max])
    title('gamma2')

    %Covariance - rho2 - at PilotPoint locations
    
    subplot(2,7,6)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_pilotpoint=1:length(V_X_ClimateDivisions)
        my_val=M_ParamsMatern_ClimateDivisions(my_type,my_pilotpoint,3);
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Y_ClimateDivisions(my_pilotpoint)/111,V_X_ClimateDivisions(my_pilotpoint)/111,'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('rho2')
    

    %Covariance - nu - at PilotPoint locations
    subplot(2,7,7)
    hold on
    my_colormap=colormap('jet');
    Plot_coastline()
    for my_pilotpoint=1:length(V_X_ClimateDivisions)
        my_max=1.5;
        my_val=M_ParamsMatern_ClimateDivisions(my_type,my_pilotpoint,4);
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Y_ClimateDivisions(my_pilotpoint)/111,V_X_ClimateDivisions(my_pilotpoint)/111,'o','MarkerSize',5,'color','k','MarkerFaceColor',my_colormap(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('nu')

    %Marginal - a0 - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamMarginal_SimulationPoints(1:length(V_X_SimulationGrid),my_type,1);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,8)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    caxis([-3 3])

    title('a0 interp')

    %Marginal - k - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamMarginal_SimulationPoints(1:length(V_X_SimulationGrid),my_type,2);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,9)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('k interp')
    caxis([0 1.5])

    %Marginal - teta - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamMarginal_SimulationPoints(1:length(V_X_SimulationGrid),my_type,3);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,10)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('teta interp')
    caxis([0 25])

    %Covariance - gamma1 - on SimulationGrid
    my_max=150;
    %my_max=40;
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamCovariance_SimulationPoints(1:length(V_X_SimulationGrid),my_type,1);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,11)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('gamma1 interp')
    caxis([-my_max my_max])

    %Covariance - gamma2 - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamCovariance_SimulationPoints(1:length(V_X_SimulationGrid),my_type,2);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,12)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('gamma2 interp')
    caxis([-my_max my_max])

    %Covariance - rho2 - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamCovariance_SimulationPoints(1:length(V_X_SimulationGrid),my_type,3);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,13)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('rho2 interp')
    caxis([0 my_max])

    %Covariance - nu - on SimulationGrid
    Mamat=[];
    ind=1;
    V_Z_SimulationGrid=M_ParamCovariance_SimulationPoints(1:length(V_X_SimulationGrid),my_type,4);
    [Mamat]=Transform_SimGrid_VectToMat(V_Z_SimulationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
    subplot(2,7,14)
    Mamat=flipud(Mamat);
    imagesc(Mamat)
    axis equal tight
    colorbar
    title('nu - interp')
    caxis([0 1.5])
end