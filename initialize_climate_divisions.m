function[nb_climate_division,M_inds_neighboors_ClimateDivisions,V_X_ClimateDivisions,V_Y_ClimateDivisions]=initialize_climate_divisions(V_Lat,V_Lon)

nb_climate_division=6;

inds_Windward_Kohala=find( (V_Lat>=19.94 & V_Lon>-155.789) | V_Lat>20.2);
inds_Leeward_Kohala=find( V_Lat<20.105 & V_Lat>19.67 & V_Lon<-155.73);
inds_Kona=find( V_Lat<19.67 & V_Lon<-155.79);
inds_Kau=find( V_Lat<19.36 & V_Lon>-155.79);
inds_Hilo=find( V_Lat>19.36 & V_Lat<19.95 & V_Lon>-155.39);
inds_Mauka=find( V_Lon>-155.78 & V_Lon<-155.39 & V_Lat<19.95 & V_Lat>19.36);

M_inds_neighboors_ClimateDivisions=NaN(6,100);
M_inds_neighboors_ClimateDivisions(1,1:length(inds_Windward_Kohala))=inds_Windward_Kohala;
M_inds_neighboors_ClimateDivisions(2,1:length(inds_Leeward_Kohala))=inds_Leeward_Kohala;
M_inds_neighboors_ClimateDivisions(3,1:length(inds_Kona))=inds_Kona;
M_inds_neighboors_ClimateDivisions(4,1:length(inds_Kau))=inds_Kau;
M_inds_neighboors_ClimateDivisions(5,1:length(inds_Hilo))=inds_Hilo;
M_inds_neighboors_ClimateDivisions(6,1:length(inds_Mauka))=inds_Mauka;

V_X_ClimateDivisions=[mean(V_Lon(inds_Windward_Kohala)); mean(V_Lon(inds_Leeward_Kohala)); mean(V_Lon(inds_Kona)); mean(V_Lon(inds_Kau)); mean(V_Lon(inds_Hilo)); mean(V_Lon(inds_Mauka)) ]*111;
V_Y_ClimateDivisions=[mean(V_Lat(inds_Windward_Kohala)); mean(V_Lat(inds_Leeward_Kohala)); mean(V_Lat(inds_Kona)); mean(V_Lat(inds_Kau)); mean(V_Lat(inds_Hilo)); mean(V_Lat(inds_Mauka)) ]*111;

end
