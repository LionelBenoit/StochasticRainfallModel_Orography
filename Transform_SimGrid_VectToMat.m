function[Mat_SimulationGrid]=Transform_SimGrid_VectToMat(V_Z_SimilationGrid,V_convert_simgrid,dimX_simGrid,dimY_simGrid,Xmin_simGrid,Xmax_simGrid,Ymin_simGrid,Ymax_simGrid,step_simulation_grid)

V_SimulationGrid=NaN(dimX_simGrid*dimY_simGrid,1);
V_SimulationGrid(V_convert_simgrid)=V_Z_SimilationGrid;
ind=1;
Mat_SimulationGrid=[];
for my_X=Xmin_simGrid:step_simulation_grid*111:Xmax_simGrid
    V_tmp=[];
    for my_Y=Ymin_simGrid:step_simulation_grid*111:Ymax_simGrid
        
        V_tmp=[V_tmp;V_SimulationGrid(ind)];
        
        ind=ind+1;
    end
    Mat_SimulationGrid=[Mat_SimulationGrid V_tmp];
end


end %end function