function[V_Z_targetpoints,alpha_TPS]=TPS_Interp(V_X_datapoints,V_Y_datapoints,V_Z_datapoints,V_X_targetpoints,V_Y_targetpoints,min_constr,max_constr,alpha_TPS,bool_plot)

if alpha_TPS==0
    %K-folds cross-validation to select alpha_thin_plate
    nb_folds=10;
    V_inds_randperm=randperm(length(V_X_datapoints));
    step_KF=floor(length(V_X_datapoints)/nb_folds);

    M_MAE_RMSE=[];

    for my_alpha=0:20:1500
        V_sim=[];
        V_obs=[];
        for my_test_set=1:nb_folds
            inds_test=V_inds_randperm(1+(my_test_set-1)*step_KF:my_test_set*step_KF);
            inds_training=V_inds_randperm;
            inds_training(inds_test)=[];

            [my_sim]=ThinPlateSmoothing(V_X_datapoints(inds_training),V_Y_datapoints(inds_training),V_Z_datapoints(inds_training),V_X_datapoints(inds_test),V_Y_targetpoints(inds_test),my_alpha);
            V_sim=[V_sim;my_sim];
            my_obs=V_Z_datapoints(inds_test);
            V_obs=[V_obs;my_obs];

        end
        MAE=mean(abs(V_sim-V_obs));
        RMSE=sqrt(mean((V_sim-V_obs).^2));

        M_MAE_RMSE=[M_MAE_RMSE;[my_alpha MAE RMSE]];
    end
    V_MAE=M_MAE_RMSE(:,2);
    ind_minMAE=find(V_MAE==min(V_MAE));
    ind_minMAE=ind_minMAE(1);
    alpha_TPS=M_MAE_RMSE(ind_minMAE,1);

    if bool_plot==1
        figure()
        subplot(1,2,1)
        plot(M_MAE_RMSE(:,1),M_MAE_RMSE(:,2),'k-')
        hold on
        plot(M_MAE_RMSE(ind_minMAE,1),M_MAE_RMSE(ind_minMAE,2),'r+')
        title('MAE=f(alphaTPS)')

        subplot(1,2,2)
        plot(M_MAE_RMSE(:,1),M_MAE_RMSE(:,3),'k-')
        title('RMSE=f(alphaTPS)')
    end
    
end
V_Z_targetpoints=ThinPlateSmoothing(V_X_datapoints,V_Y_datapoints,V_Z_datapoints,V_X_targetpoints,V_Y_targetpoints,alpha_TPS);
V_Z_targetpoints=max(min(V_Z_targetpoints,max_constr),min_constr);
end