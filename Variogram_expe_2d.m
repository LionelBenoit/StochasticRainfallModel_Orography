function [meanh,gammah,nbpairs]=Variogram_expe_2d(pts,nblags,maxh,sampling)

% See DS user guide for functions documentation

% Written by Gregoire Mariethoz, 2010 - modified LB 2022
xy=pts(:,[1 2]);
v=pts(:,3);

ind=randperm(size(xy,1));
xy=xy(ind(1:sampling:end),:);
v=v(ind(1:sampling:end));

meanh=zeros(1,nblags);
gammah=zeros(1,nblags);

%calculating distances matrix
h=squareform(pdist(xy));
gamma=0.5*squareform(pdist(v).^2);

step=maxh/nblags;
%looping categories
for i=1:nblags
    maxlag=i*step;
    category_points=(h<=maxlag & h>(i-1)*step);
    meanh(i)=mean(mean(h(category_points)));
    gammah(i)=mean(gamma(category_points));
    nbpairs(i)=sum(category_points(:));
end

