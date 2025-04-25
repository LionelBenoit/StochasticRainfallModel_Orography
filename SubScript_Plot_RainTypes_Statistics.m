V_Raintypes_Plot=V_RainTypes_Training;
M_data_Plot=M_Data_Training;
M_datevec_Plot=M_datevec_Training;

my_markersize=8;
for my_type=0:nb_clusters

    figure(num_fig+my_type)
    clf

    if my_type==0
        inds_type=find(V_Raintypes_Plot>0);
    else
        inds_type=find(V_Raintypes_Plot==my_type);
    end

    V_monthly_freq=zeros(12,1);
    for my_month=1:12
        V_monthly_freq(my_month)=sum(M_datevec_Plot(inds_type,2)==my_month)/sum(M_datevec_Plot(:,2)==my_month)*100;
    end
    subplot(1,6,1)
    plot(1:12,V_monthly_freq,'k-','LineWidth',2)
    axis([0.5 12.5 0 70])
    title('Monthly freq. of occurence (%)')
    
    subplot(1,6,2)
    my_colormapjet=colormap('jet');
    Plot_coastline()

    my_max=15;
    for my_gauge=1:length(V_Lon)
        my_val=mean(M_data_Plot(my_gauge,inds_type),'omitnan');
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',my_markersize,'color','k','MarkerFaceColor',my_colormapjet(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('E(R) [mm/day]')
    %---
    subplot(1,6,3)
    my_colormapjet=colormap('jet');
    Plot_coastline()

    my_max=100;
    for my_gauge=1:length(V_Lon)
        my_val=sum(M_data_Plot(my_gauge,inds_type)>0)/length(inds_type)*100;
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',my_markersize,'color','k','MarkerFaceColor',my_colormapjet(my_ind_color,:))
    end
    colorbar
    caxis([0 my_max])
    title('P(R>0) [%]')

    subplot(1,6,4)
    my_colormapjet=colormap('jet');
    Plot_coastline()

    my_max=4;
    for my_gauge=1:length(V_Lon)
        inds_select=find(M_data_Plot(my_gauge,inds_type)>0.1);
        my_val=mean(M_data_Plot(my_gauge,inds_type(inds_select)));
        my_val=log(my_val+1);
        my_ind_color=min(floor(my_val/my_max*255)+1,255);
        plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',my_markersize,'color','k','MarkerFaceColor',my_colormapjet(my_ind_color,:))
    end

    cbh=colorbar;
    cbh.Ticks=[0 log(3) log(8) log(21) log(51)];
    cbh.TickLabels={'0','2','7','20','50'};
    hold on
    caxis([0 my_max])
    title('E(R|R>0) [mm/day]')


    subplot(1,6,5)
    my_colormapjet=colormap('jet');
    Plot_coastline()
    
    Lat_target=19.55;
    Lon_target=-155.91;
    V_dist=sqrt( (V_Lon-Lon_target).^2 + (V_Lat-Lat_target).^2);
    ind_target_gauge=find(V_dist == min(V_dist));

    my_max=1;
    for my_gauge=1:length(V_Lon)
        if my_gauge==ind_target_gauge
            plotm(V_Lat(my_gauge),V_Lon(my_gauge),'*','MarkerSize',my_markersize,'color','k')
        else
            inds_corr=find(M_data_Plot(my_gauge,inds_type)>0.1 & M_data_Plot(ind_target_gauge,inds_type)>0.1);
            if length(inds_corr)>30
                my_val=corr(M_data_Plot(my_gauge,inds_type(inds_corr))',M_data_Plot(ind_target_gauge,inds_type(inds_corr))','Type','Pearson');
                my_val=max(my_val,0);
                my_ind_color=min(floor(my_val/my_max*255)+1,255);
                plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',my_markersize,'color','k','MarkerFaceColor',my_colormapjet(my_ind_color,:))
            else
                plotm(V_Lat(my_gauge),V_Lon(my_gauge),'.','MarkerSize',my_markersize,'color','k')
            end
        end
    end
    plotm(V_Lat(ind_target_gauge),V_Lon(ind_target_gauge),'*','MarkerSize',12,'color','k')
    
    cbh=colorbar;
    cbh.Ticks=[0 0.5 1];
    cbh.TickLabels={'0','0.5','1'};
    caxis([0.0 my_max])
    title('Pearson corr(R>0) with *')


    subplot(1,6,6)
    my_colormapjet=colormap('jet');
    Plot_coastline()
    
    Lat_target=19.60;
    Lon_target=-155.08;
    V_dist=sqrt( (V_Lon-Lon_target).^2 + (V_Lat-Lat_target).^2);
    ind_target_gauge=find(V_dist == min(V_dist));

    my_max=1;
    for my_gauge=1:length(V_Lon)
        if my_gauge==ind_target_gauge
            plotm(V_Lat(my_gauge),V_Lon(my_gauge),'*','MarkerSize',my_markersize,'color','k')
        else
            inds_corr=find(M_data_Plot(my_gauge,inds_type)>0.1 & M_data_Plot(ind_target_gauge,inds_type)>0.1);
            if length(inds_corr)>30
                my_val=corr(M_data_Plot(my_gauge,inds_type(inds_corr))',M_data_Plot(ind_target_gauge,inds_type(inds_corr))','Type','Pearson');
                my_val=max(my_val,0);
                my_ind_color=min(floor(my_val/my_max*255)+1,255);
                plotm(V_Lat(my_gauge),V_Lon(my_gauge),'o','MarkerSize',my_markersize,'color','k','MarkerFaceColor',my_colormapjet(my_ind_color,:))
            else
                plotm(V_Lat(my_gauge),V_Lon(my_gauge),'.','MarkerSize',my_markersize,'color','k')
            end
        end
    end
    plotm(V_Lat(ind_target_gauge),V_Lon(ind_target_gauge),'*','MarkerSize',12,'color','k')

    
    cbh=colorbar;
    cbh.Ticks=[0 0.5 1];
    cbh.TickLabels={'0','0.5','1'};
    caxis([0 my_max])
    title('Pearson corr(R>0) with *')
   
end