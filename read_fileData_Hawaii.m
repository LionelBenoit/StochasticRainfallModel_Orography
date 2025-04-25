function[V_Time,M_Data_TargetIsland,V_Lon_TargetIsland,V_Lat_TargetIsland,V_Alti_TargetIsland]=read_fileData_Hawaii(my_filename)

M = readtable(my_filename,'Delimiter',',','FileType','text');


V_Lat=M.LAT;
V_Lon=M.LON;
V_alti=M.ELEV_m_;

V_Inds_BigIsland=[];

M_Data=NaN(size(M,1),size(M,2)-14);

V_Island=M(:,5);
for i=1:size(M_Data,1)
    if strcmp(V_Island{i,V_Island.Properties.VariableNames},'BI')==1
        V_Inds_BigIsland=[V_Inds_BigIsland; i];
    end
end

for my_col=14:size(M,2)
    my_col
    my_name=M.Properties.VariableNames{my_col};
    my_data_col=M.(my_name);
    if strcmp(class(my_data_col),'cell')==1
        for i=1:length(my_data_col)
            if strcmp(my_data_col{i},"NA")==0
                M_Data(i,my_col-14+1)=str2double(my_data_col{i});
            end
        end
    else
        M_Data(:,my_col-14+1)=my_data_col;
        
    end
end

M_Data_BigIsland=M_Data(V_Inds_BigIsland,:);


%Select area of interest and remove problematic gauges
V_Lon_BigIsland=V_Lon(V_Inds_BigIsland);
V_Lat_BigIsland=V_Lat(V_Inds_BigIsland);
V_Alti_BigIsland=V_alti(V_Inds_BigIsland);

V_Lon_TargetIsland=V_Lon_BigIsland;
V_Lat_TargetIsland=V_Lat_BigIsland;
V_Alti_TargetIsland=V_Alti_BigIsland;

M_Data_TargetIsland=M_Data_BigIsland(:,366:end);%#ToSet

t0=datenum('01-Jan-1991','dd-mmm-yyyy');%#ToSet
V_Time=t0:t0+size(M_Data_TargetIsland,2);
V_Time=V_Time';

M_Data_TargetIsland=M_Data_TargetIsland(:,3288:10592);
V_Time=V_Time(3288:10592);
end