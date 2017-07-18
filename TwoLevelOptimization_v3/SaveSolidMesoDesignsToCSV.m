function []=SaveSolidMesoDesignsToCSV() 

ne=50*20;
folderNum=0;
macro_meso_iteration=1;

mesoDesign = ones(35,35);

for i=1:ne
    elementNumber=i;
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
%     outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,settings.macro_meso_iteration,elementNumber);
    x = mesoDesign;
    csvwrite(outname,x);
end
