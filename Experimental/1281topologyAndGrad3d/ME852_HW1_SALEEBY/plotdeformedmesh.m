clear all;
clf;

load nodes.dat;
load elements.dat;
load feU.dat;

feU=feU*80;

nNodes=size(nodes,1);
nodes=nodes(:,2:4);

for i=1:1:nNodes;
    nodes(i,1)=nodes(i,1)+feU(i,5);
    nodes(i,2)=nodes(i,2)+feU(i,6);
    nodes(i,3)=nodes(i,3)+feU(i,7);
end 

p=patch('Vertices',nodes,'Faces',elements(:,2:5));
set(p,'facecolor','green','edgecolor','black');

s1=elements(:,[2 3 7 6]);
p=patch('Vertices',nodes,'Faces',s1);
set(p,'facecolor','green','edgecolor','black');


s1=elements(:,[3 4 8 7]);
p=patch('Vertices',nodes,'Faces',s1);
set(p,'facecolor','green','edgecolor','black');


s1=elements(:,[4 5 9 8]);
p=patch('Vertices',nodes,'Faces',s1);
set(p,'facecolor','green','edgecolor','black');


s1=elements(:,[5 2 6 9]);
p=patch('Vertices',nodes,'Faces',s1);
set(p,'facecolor','green','edgecolor','black');

p=patch('Vertices',nodes,'Faces',elements(:,6:9));
set(p,'facecolor','green','edgecolor','black');

daspect([1 1 0.7]);
view(30,50); 
grid on;
camlight; lighting gouraud;

alpha(.75)

axis([-60 60 -10 70 -5 15]);


print -r200 -djpeg90 deformedmesh.jpg


