%%---- VJC'S VERSION OF "bwdist" ----
function dist = bwdist(image)
% Gives Euclidean distance to the closest non-zero entry in "image"
dist=zeros(size(image));
[yinds,xinds]=find(image~=0);
for i=1:size(image,2)
 for j=1:size(image,1)
  dist(j,i)=sqrt(min((yinds-j).^2+(xinds-i).^2));
 end
end
