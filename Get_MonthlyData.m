function[V_X_MonthlyCondi,V_Y_MonthlyCondi,V_Obs_MonthlyCondi]=Get_MonthlyData(my_year,my_month,V_X_ObsSites,V_Y_ObsSites,V_X_SimulationGrid,V_Y_SimulationGrid,nb_condiMonthly)

[V_X_MonthlyCondi,V_Y_MonthlyCondi]=Set_ConditioningLocations_Monthly(V_X_ObsSites,V_Y_ObsSites,V_X_SimulationGrid,V_Y_SimulationGrid,nb_condiMonthly);

%set target month and get the monthly map

cd('Data')
cd('BI_MONTHLY_2018-2020')
my_file=strcat(num2str(my_year),'_',num2str(my_month,'%02d'),'_bi_rf_mm.tif');
[A,R]=readgeoraster(my_file);
cd('..')
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

%Extract monthly data at pseudo-observation locations
V_Obs_MonthlyCondi=NaN(size(V_X_MonthlyCondi));

for i=1:length(V_X_MonthlyCondi)
    V_dist= sqrt( (V_X_MonthlyCondi(i)-V_X_Monthly_HR).^2 + (V_Y_MonthlyCondi(i)-V_Y_Monthly_HR).^2 );
    dist_min= min(V_dist);
    ind=find(V_dist==dist_min);
    ind=ind(1);
    V_Obs_MonthlyCondi(i)=V_Obs_Monthly_HR(ind);
    V_X_MonthlyCondi(i)=V_X_Monthly_HR(ind);
    V_Y_MonthlyCondi(i)=V_Y_Monthly_HR(ind);
end

%---
V_DistToSimGrid=NaN(size(V_X_MonthlyCondi));
for i=1:length(V_DistToSimGrid)
    V_dist=sqrt( (V_X_MonthlyCondi(i)-V_X_SimulationGrid).^2 + (V_Y_MonthlyCondi(i)-V_Y_SimulationGrid).^2 );
    V_DistToSimGrid(i)=min(V_dist);
end

V_X_MonthlyCondi(V_DistToSimGrid<0.01)=[];
V_Y_MonthlyCondi(V_DistToSimGrid<0.01)=[];
V_Obs_MonthlyCondi(V_DistToSimGrid<0.01)=[];

end