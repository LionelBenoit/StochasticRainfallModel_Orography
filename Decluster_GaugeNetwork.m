function[inds_decluster]=Decluster_GaugeNetwork(dist_decluster,V_Lon,V_Lat,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid)

%Keep 1 gauge by gridcell of size dist_decluster x dist_decluster
%Keep the closest gauge to the center of the gridcell

% dist_decluster is in km

inds_decluster=[];
for my_X=Xmin_simGrid-dist_decluster/2:dist_decluster:Xmax_simGrid+dist_decluster/2
    for my_Y=Ymin_simGrid-dist_decluster/2:dist_decluster:Ymax_simGrid+dist_decluster/2
        V_dist=sqrt((V_Lon*111-my_X).^2+(V_Lat*111-my_Y).^2);
        dist_min=min(V_dist);
        if dist_min<dist_decluster/2
            ind_closest=find(V_dist==dist_min);
            inds_decluster=[inds_decluster;ind_closest];
        end
    end
end
inds_decluster=unique(inds_decluster);

end %end function