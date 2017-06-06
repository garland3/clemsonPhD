function [DV] = GenerateMesoValidationTargets(DV,config,matProp)
fprintf('Generate target values for meso validation mode\n');
Usize = size(config.loadingCase,2);
nnodes =(config.nelx+1)*(config.nely+1)*2;
DV.U=ones(Usize,nnodes);% 1; % so it is not NULL. Just assign a placeholder value
DV.sensitivityElastic=1;
DV.sensitivityElasticPart2=1;
DV.w=1;
DV.lambda1=1;

% DV.x = ones(config.nely, config.nelx);

% totalValidationProblems=config.nely*config.nelx;
% numSegments = floor(totalValidationProblems^(1/3));
% numSegmentsExx = numSegments;
% numSegmentsTheta = floor(totalValidationProblems/(numSegmentsExx^2));
% 
% numSegmentsExx=numSegmentsExx-1;
% % numSegmentsTheta=numSegmentsTheta-1;

temp=config.validationGridSizeNelx-1;

ExxVector =0:matProp.E_material1/temp:matProp.E_material1;
EyyVector =0:matProp.E_material1/temp:matProp.E_material1;
thetaVector = 0:(pi/2)/temp:pi/2;
[ExxValues, EyyValues, ThetaValues] = meshgrid(ExxVector,EyyVector,thetaVector);

[ t1 t2 t3]=size(ExxValues);
ne = t1*t2*t3;
simpDensity = ones(ne,1);

% xValues
% % save the density field
% outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
% csvwrite(outname,DV.x);
folderNum=0;
% save the Exx field
mm_iteration=1;
outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,ExxValues);

% save the Eyy field
outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,EyyValues);

% save the Theta field
outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,ThetaValues);

% Save a placeholder file. 
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
csvwrite(outname,[1 2]);

outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
csvwrite(outname,simpDensity);

%ExxValues=padarray(ExxValues,
% ExxValues=reshape(ExxValues,1,[]);
% [t1 t2]=size(ExxValues);
% ExxValues=padarray(ExxValues,[0 totalValidationProblems-t2],'post');
% ExxValues=reshape(ExxValues,config.nely,config.nelx);
% 
% EyyValues=reshape(EyyValues,1,[]);
% EyyValues=padarray(EyyValues,[0 totalValidationProblems-t2],'post');
% EyyValues=reshape(EyyValues,config.nely,config.nelx);
% 
% ThetaValues=reshape(ThetaValues,1,[]);
% ThetaValues=padarray(ThetaValues,[0 totalValidationProblems-t2],'post');
% ThetaValues=reshape(ThetaValues,config.nely,config.nelx);
% 
% xSimp = ones(t1, t2);
% xSimp=reshape(xSimp,1,[]);
% xSimp=padarray(xSimp,[0 totalValidationProblems-t2],'post');
% xSimp=reshape(xSimp,config.nely,config.nelx);
% 
% DV.Exx = ExxValues;
% DV.Eyy = EyyValues;
% DV.t = ThetaValues;
% DV.x = xSimp;

% SaveMacroProblemStateToCSV(config,DV,matProp)

end
