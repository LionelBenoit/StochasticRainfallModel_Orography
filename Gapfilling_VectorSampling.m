function[M_Data2,V_Lon2,V_Lat2]=Gapfilling_VectorSampling(M_Data,V_Lon,V_Lat,threshold_to_be_filled)

M_Data_TargetIsland=M_Data;

M_Isnan=isnan(M_Data_TargetIsland);
V_workingRatio_Stations=1-sum(M_Isnan,2)/size(M_Isnan,2);

%{
figure(1)
clf
imagesc(M_Isnan)
colorbar

figure(2)
clf
imagesc(V_workingRatio_Stations)
colorbar
%}
V_ind_gauges_toselect=(1-sum(M_Isnan,2)/size(M_Isnan,2))>threshold_to_be_filled;%#ToSet


V_Lon2=V_Lon(V_ind_gauges_toselect);
V_Lat2=V_Lat(V_ind_gauges_toselect);
M_Data2=M_Data(V_ind_gauges_toselect,:);

%gap filling
M3=M_Data2';
M3_filled=M_Data2';
for my_station=1:size(M3,2)
    my_station
    my_inds_to_fill=find(isnan(M3(:,my_station)));
    M_guiding=[M3(~isnan(M3(:,my_station)),1:my_station-1) M3(~isnan(M3(:,my_station)),my_station+1:end)];
    V_select=M3(~isnan(M3(:,my_station)),my_station);
    for j=1:length(my_inds_to_fill)
        Mamat=M_guiding;
        V_target=[M3(my_inds_to_fill(j),1:my_station-1) M3(my_inds_to_fill(j),my_station+1:end)];
        
        V_dist=nanmean(abs(Mamat-repmat(V_target,size(Mamat,1),1)),2);
        [~,ind_to_pick]=min(V_dist);
        ind_to_pick=ind_to_pick(1);
        M3_filled(my_inds_to_fill(j),my_station)=V_select(ind_to_pick);
    end
end

M_Data2=M3_filled';
end