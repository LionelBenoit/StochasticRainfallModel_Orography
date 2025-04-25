function[V_X_SimilationGrid,V_Y_SimilationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid]=Create_SimulationGrid_fromGeoRaster(A,R,step_simulation_grid)

V_X_SimilationGrid=[];
V_Y_SimilationGrid=[];
V_convert_simgrid=[];
my_ind=1;

for my_X=min(R.LongitudeLimits)*111:step_simulation_grid*111:max(R.LongitudeLimits)*111
    for my_Y=min(R.LatitudeLimits)*111:step_simulation_grid*111:max(R.LatitudeLimits)*111
        

        my_l=min(floor((max(R.LatitudeLimits)*111-my_Y)/(max(R.LatitudeLimits)*111 - min(R.LatitudeLimits)*111)*size(A,1))+1 , size(A,1));
        my_c=min(floor((my_X-min(R.LongitudeLimits)*111)/(max(R.LongitudeLimits)*111 - min(R.LongitudeLimits)*111)*size(A,2))+1 , size(A,2));

        if A(my_l,my_c)>-999
            V_X_SimilationGrid=[V_X_SimilationGrid;my_X];
            V_Y_SimilationGrid=[V_Y_SimilationGrid;my_Y];
            V_convert_simgrid=[V_convert_simgrid;my_ind];
        end

        my_ind=my_ind+1;

    end
end
dimX_simGrid=length(min(R.LongitudeLimits)*111:step_simulation_grid*111:max(R.LongitudeLimits)*111);
dimY_simGrid=length(min(R.LatitudeLimits)*111:step_simulation_grid*111:max(R.LatitudeLimits)*111);
Xmin_simGrid=min(R.LongitudeLimits)*111;
Xmax_simGrid=max(R.LongitudeLimits)*111;
Ymin_simGrid=min(R.LatitudeLimits)*111;
Ymax_simGrid=max(R.LatitudeLimits)*111;
end %end function