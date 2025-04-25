figure(100)
clf

Plot_coastline()
hold on
for my_VirtualStation=1:length(V_Y_MonthlyCondi)
    plotm(V_Y_MonthlyCondi(my_VirtualStation)/111,V_X_MonthlyCondi(my_VirtualStation)/111,'+','MarkerSize',7,'color','k')
    
end
for my_gauge=1:length(V_Lat)
    plotm(V_Lat(my_gauge),V_Lon(my_gauge),'*','MarkerSize',7,'color','k')
   
end
title('Network of rain gauges (*) and virtual stations (+)')

figure(101)
clf

V_X_Condi=V_Lon*111;
V_Y_Condi=V_Lat*111;

X_min=Xmin_simGrid;
X_max=Xmax_simGrid;
Y_min=Ymin_simGrid;
Y_max=Ymax_simGrid;

my_markersize=5;
my_colormap=colormap('jet');
my_colormap(256,:)=[1 1 1];

V_condi_data=M_Data_Month(:,1);
V_monthlySum_Obs=zeros(size(V_condi_data));

for my_day=1:length(inds_days_month)
    V_condi_data=M_Data_Month(:,my_day); 
    V_monthlySum_Obs=V_monthlySum_Obs+V_condi_data;
end

my_max=950;

% Subplot(1,3,1): Load and plot monthly rain map from HCDP
colormap('jet');
cd('BI_MONTHLY_2018-2020')
my_file=strcat(num2str(my_year),'_',num2str(my_month,'%02d'),'_bi_rf_mm.tif');
[A,R]=readgeoraster(my_file);
cd('..')

V_Lon_Monthly_HR=NaN(size(A,1)*size(A,2),1);
V_Lat_Monthly_HR=NaN(size(A,1)*size(A,2),1);
V_Obs_Monthly_HR=NaN(size(A,1)*size(A,2),1);
ind=1;

for my_mat_line=1:size(A,1)
    for my_mat_col=1:size(A,2)
        if A(my_mat_line,my_mat_col)>=0
            V_Obs_Monthly_HR(ind)=A(my_mat_line,my_mat_col);
            V_Lat_Monthly_HR(ind)=max(R.LatitudeLimits)-(my_mat_line-1)*R.CellExtentInLatitude;
            V_Lon_Monthly_HR(ind)=min(R.LongitudeLimits)+(my_mat_col-1)*R.CellExtentInLongitude;
            ind=ind+1;
        end
    end
end
V_Obs_Monthly_HR=V_Obs_Monthly_HR(1:ind-1);
V_Lon_Monthly_HR=V_Lon_Monthly_HR(1:ind-1);
V_Lat_Monthly_HR=V_Lat_Monthly_HR(1:ind-1);

V_X_Monthly_HR=V_Lon_Monthly_HR*111;
V_Y_Monthly_HR=V_Lat_Monthly_HR*111;

V_Obs_MonthlyHCDP=NaN(size(V_X_SimulationGrid));

for i=1:length(V_X_SimulationGrid)
    V_dist= sqrt( (V_X_SimulationGrid(i)-V_X_Monthly_HR).^2 + (V_Y_SimulationGrid(i)-V_Y_Monthly_HR).^2 );
    dist_min= min(V_dist);
    ind=find(V_dist==dist_min);
    ind=ind(1);
    V_Obs_MonthlyHCDP(i)=min(V_Obs_Monthly_HR(ind),my_max-10);
end

V_valueToPlot=V_Obs_MonthlyHCDP;
[Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);

Mamat(isnan(Mamat))=my_max+1;

colormap('jet');
subplot(1,3,1)

imagesc(flipud(Mamat))
axis equal tight
clim([0 my_max])
colorbar
title('Reference Monthly Rainfall')



% Subplot(1,3,2): plot monthly rain map from condition simulation - conditioning to daily values only

%initialization
M_simul_condi=M_SimulCondi_DailyOnly(:,:,1);
[Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
M_monthlySum_Map=zeros(size(Mamat));

for my_day=1:length(inds_days_month)
    if V_RainTypes_Month(my_day)>0
        M_simul_condi=M_SimulCondi_DailyOnly(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        [Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
        M_monthlySum_Map=M_monthlySum_Map+Mamat;
    end
end
M_monthlySum_Map(isnan(M_monthlySum_Map))=my_max+1;

colormap('jet');
subplot(1,3,2)
imagesc(flipud(M_monthlySum_Map))
axis equal tight
clim([0 my_max])
colorbar
title('Monthly Rainfall - gauges only')

%----
% Subplot(1,3,3): plot monthly rain map from condition simulation - conditioning to daily values and monthly totals at selected locations

%initialization
M_simul_condi=M_SimulCondi_MwG(:,:,1);
[Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
M_monthlySum_Map=zeros(size(Mamat));

for my_day=1:length(inds_days_month)
    if V_RainTypes_Month(my_day)>0
        M_simul_condi=M_SimulCondi_MwG(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        [Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
        M_monthlySum_Map=M_monthlySum_Map+Mamat;
    end
end
M_monthlySum_Map(isnan(M_monthlySum_Map))=my_max+1;

colormap('jet');
subplot(1,3,3)
imagesc(flipud(M_monthlySum_Map))
axis equal tight
clim([0 my_max])
colorbar
title('Monthly Rainfall - gauges+virtual stations')

colormap('jet');
%---
hold on
for my_gauge=1:length(V_X_Condi)
    my_X=(V_X_Condi(my_gauge)-X_min+step_simulation_grid*111)/(step_simulation_grid*111);
    my_Y=(Y_max-V_Y_Condi(my_gauge)+step_simulation_grid*111)/(step_simulation_grid*111);
    plot(my_X,my_Y,'o','MarkerSize',my_markersize,'color','w')
end

for my_VirtualStation=1:length(V_X_Condi)
    my_X=(V_X_MonthlyCondi(my_VirtualStation)-X_min+step_simulation_grid*111)/(step_simulation_grid*111);
    my_Y=(Y_max-V_Y_MonthlyCondi(my_VirtualStation)+step_simulation_grid*111)/(step_simulation_grid*111);
    plot(my_X,my_Y,'+','MarkerSize',my_markersize,'color','w')
end

%}

%%
%Plot daily rain maps

for my_day=1:length(inds_days_month)
    f=figure(1000+my_day);
    clf
    my_colormap=colormap('jet');
    my_colormap(256,:)=[1 1 1];
    f.Position = [100 100 1000 400];

    V_condi_data=M_Data_Month(:,my_day); 

    if V_RainTypes_Month(my_day)>0
        my_max1=max(V_condi_data);

        M_simul_condi=M_SimulCondi_DailyOnly(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        my_max2=max(V_valueToPlot);

        M_simul_condi=M_SimulCondi_MwG(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        my_max3=max(V_valueToPlot);

        my_max=max([my_max1 my_max2 my_max3]);

        ax=subplot(1,2,1);
        M_simul_condi=M_SimulCondi_DailyOnly(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        [Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
        Mamat(isnan(Mamat))=my_max+1;

        dim = [0.275 0.02 0.3 0.3];
        str = {'Rain gauges only'};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
        colormap(my_colormap)
        imagesc(flipud(Mamat))
        axis equal tight
        clim([0 my_max+1])
        colorbar
        title(strcat(num2str(M_datevec_Test(inds_days_month(my_day),1)),'/',num2str(M_datevec_Test(inds_days_month(my_day),2),'%02.0f'),'/',num2str(M_datevec_Test(inds_days_month(my_day),3),'%02.0f'),' - RainType:',num2str(V_RainTypes_Test(inds_days_month(my_day)))))
        
        hold on
        for my_gauge=1:length(V_X_Condi)
            my_val=V_condi_data(my_gauge);
            my_X=(V_X_Condi(my_gauge)-X_min+step_simulation_grid*111)/(step_simulation_grid*111);
            my_Y=(Y_max-V_Y_Condi(my_gauge)+step_simulation_grid*111)/(step_simulation_grid*111);
            my_ind_color=min(floor(my_val/my_max*255)+1,255);
            plot(my_X,my_Y,'o','MarkerSize',my_markersize,'color','w','MarkerFaceColor',my_colormap(my_ind_color,:))
        end
        set(ax,'xticklabel',[])
        set(ax,'yticklabel',[])

        ax=subplot(1,2,2);
        M_simul_condi=M_SimulCondi_MwG(:,:,my_day);
        V_valueToPlot=median(M_simul_condi,2);
        [Mamat]=Transform_SimGrid_VectToMat(V_valueToPlot,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid);
        Mamat(isnan(Mamat))=my_max+1;
        
        colormap(my_colormap)
        imagesc(flipud(Mamat))
        axis equal tight
        clim([0 my_max+1])
        colorbar
        
        hold on
        for my_gauge=1:length(V_X_Condi)
            my_val=V_condi_data(my_gauge);
            my_X=(V_X_Condi(my_gauge)-X_min+step_simulation_grid*111)/(step_simulation_grid*111);
            my_Y=(Y_max-V_Y_Condi(my_gauge)+step_simulation_grid*111)/(step_simulation_grid*111);
            my_ind_color=min(floor(my_val/my_max*255)+1,255);
            plot(my_X,my_Y,'o','MarkerSize',my_markersize,'color','w','MarkerFaceColor',my_colormap(my_ind_color,:))
        end

        dim = [0.715 0.02 0.3 0.3];
        str = {'Rain gauges','+Virtual stations'};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');

        set(ax,'xticklabel',[])
        set(ax,'yticklabel',[])
        
        %saveas(f, strcat(num2str(M_datevec_Test(inds_days_month(my_day),1)),num2str(M_datevec_Test(inds_days_month(my_day),2), '%02.0f'),num2str(M_datevec_Test(inds_days_month(my_day),3),'%02.0f'),'.png'), 'png')
    end
end
