function[]=Plot_coastline()

worldmap([18.80 20.35],[-156.20 -154.70]);
geoshow('Coastline_Hawaii.shp','FaceColor','white')
mlabel off; plabel off; gridm off

end