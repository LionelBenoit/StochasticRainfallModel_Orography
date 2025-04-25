function[V_RainTypes,my_GMModel]=RainTyping(V_p0,V_k,V_teta,M_PCs,nb_clusters)

V_RainTypes=zeros(size(V_p0));

inds_dry=(V_k==-1);
inds_wet=(V_k>0);

V_p0(inds_dry)=[];
V_k(inds_dry)=[];
V_teta(inds_dry)=[];
M_PCs(inds_dry,:)=[];

V_p0(V_p0<0.01)=0.01;
V_p0(V_p0>0.99)=0.99;

V_k(V_k<0.01)=0.01;

V_teta(V_teta<0.01)=0.01;

M_PCs=M_PCs(:,1:3);%only 3 first PCs


M_rainpatterns=[norminv(V_p0);log(V_k);log(V_teta);M_PCs'];


display('before cluster')
options = statset('MaxIter',1e5,'TolFun',1e-20);
my_GMModel = fitgmdist(M_rainpatterns',nb_clusters,'Options',options,'CovarianceType','diagonal','RegularizationValue',0.1);
display('after cluster')
[V_RainTypes_tmp,~,~] = cluster(my_GMModel,M_rainpatterns');
V_RainTypes(inds_wet)=V_RainTypes_tmp;

end
