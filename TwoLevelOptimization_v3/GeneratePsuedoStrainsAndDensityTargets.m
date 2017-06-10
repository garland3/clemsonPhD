function [] =   GeneratePsuedoStrainsAndDensityTargets(DV,config,matProp, step)

% Usize = size(config.loadingCase,2);
% nnodes=config.targetTestVectorLen^4*2;
% DV.U=ones(Usize,nnodes);% 1; % so it is not NULL. Just assign a placeholder value
% DV.sensitivityElastic=1;
% DV.sensitivityElasticPart2=1;
% DV.w=1;
% DV.lambda1=1;
temp=config.targetTestVectorLen-1;
term1 =0:1/temp:1; % -1 to 1 is the domain
term2 =0:1/temp:1;  % -1 to 1 is the domain
term3 =-0.5:2/temp:0.5;  % -1 to 1 is the domain
densityTargetsVector = 2/temp:1/temp:1-1/temp;  % 0 to 1 is the domain
[strain1, strain2, strain3, densityTargets] = ndgrid(term1,term2,term3,densityTargetsVector);
SIMPTemp = reshape(strain1,[],1);
[ numProblems t1]=size(SIMPTemp)

folderNum=0;
if(step==1)
    fprintf('Generate target values for meso validation mode\n');
    mm_iteration=1;
    outname = sprintf('./out%i/strain1%i.csv',folderNum,mm_iteration);
    csvwrite(outname,strain1);


    outname = sprintf('./out%i/strain2%i.csv',folderNum,mm_iteration);
    csvwrite(outname,strain2);

    outname = sprintf('./out%i/strain3%i.csv',folderNum,mm_iteration);
    csvwrite(outname,strain3);


    outname = sprintf('./out%i/strain3%i.csv',folderNum,mm_iteration);
    csvwrite(outname,strain3);


    outname = sprintf('./out%i/densityTargets%i.csv',folderNum,mm_iteration);
    csvwrite(outname,densityTargets);


    % Save a placeholder file. 
    outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
    csvwrite(outname,[1 2]);

    % The C program reads this array to figure out how many jobs there are. 
    
    outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
    csvwrite(outname,SIMPTemp);

elseif(step==2)
    %  ----------------------------------
    %
    %           Part 2
    % 
    %            Read the results. Map D_h 1,1
    %
    %  ----------------------------------
  
    macro_meso_iteration=1;
    ne=numProblems;
    
    folderNum=0;
    mm_iteration=1;
    % pseudoStrain
    outname = sprintf('./out%i/strain1%i.csv',folderNum,mm_iteration);
    strain1FromCSV= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain2%i.csv',folderNum,mm_iteration);
    strain2FromCSV= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain3%i.csv',folderNum,mm_iteration);
    strain3FromCSV=csvread(outname);
    
    % targetDensity
    outname = sprintf('./out%i/densityTargets%i.csv',folderNum,mm_iteration);
    densityTargetFromCSV=csvread(outname);
    
    % The when saving the 4D array, it is converted into 2D. This
    % conversion might changing the ordering of the elements, so to make
    % sure they are correct, just reread the .csv files. 
    p = plotResults;
    
    
    D11 = [];
    D12=[];
    D22=[];
    D33=[];
    
    pstrain1=[];
    pstrain2=[];
    pstrain3=[];
    etaTarget = [];
    for e = 1: ne
        fprintf('element %i of %i\n',e,ne);
        
        outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
        
        if exist(outname, 'file') ~= 2
            fprintf('File does not exist. Retry\n');
            combinedTopologyOptimization('1', '1', '1','100', int2str(e));
            %continue;
        end
        %          outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,macro_meso_iteration,elementNumber);
        try
            Din = csvread(outname);
        catch
            str = sprintf('error reading a file\n'); display(str);
            continue
        end
        
        D11Local = Din(1,1);
        D12Local = Din(1,2);
        D22Local = Din(2,2);
        D33Local = Din(3,3);
        P1_local = strain1FromCSV(e);
        p2_local =strain2FromCSV(e);
        p3_local = strain3FromCSV(e);
        etaLocal = densityTargetFromCSV(e);
        
        
        %   fprintf('D11 \t\t%f, \nD22 \t\t%f\nD33 \t\t%f\npstrain \t%f %f %f\nDensity \t%f\n\n',D11Local,D11Local,D33Local,P1_local,p2_local,p3_local,etaLocal);
        
        outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,e);
        etaArray = csvread(outname);
        
        if 1==0
            p.PlotArrayGeneric(etaArray,'eta')
            title({sprintf('pstrain %f %f %f',P1_local,p2_local,p3_local),sprintf('Density %f',etaLocal),sprintf('D11 %f, D22 %f,D33 %f',D11Local,D11Local,D33Local )  })
            caxis([0 1])
            drawnow
            
        end
        
        
        
        D11=[D11 D11Local];
        D12=[D12 D12Local];
        D22=[D22 D22Local];
        D33=[D33 D33Local];
        
        pstrain1=[pstrain1 P1_local];
        pstrain2=[pstrain2 p2_local];
        pstrain3=[pstrain3 p3_local];
        etaTarget=[etaTarget etaLocal];
        
        
    end
    
    % ---------------
    % Add a zero case 
    %-----------------   
    D11=[D11 0];
    D12=[D12 0];
    D22=[D22 0];
    D33=[D33 0];
    
    pstrain1=[pstrain1 0];
    pstrain2=[pstrain2 0];
    pstrain3=[pstrain3 0];
    etaTarget=[etaTarget 0];
    
     % ---------------
    % Add a one case case , totally solid
    %-----------------   
    D11=[D11 matProp.E_material1/(1-matProp.v^2)];
    D12=[D12 matProp.v*matProp.E_material1/(1-matProp.v^2)];
    D22=[D22 matProp.E_material1/(1-matProp.v^2)];
    D33=[D33 1/2*(1-matProp.v)*matProp.E_material1/(1-matProp.v^2)];
    
    pstrain1=[pstrain1 0];
    pstrain2=[pstrain2 0];
    pstrain3=[pstrain3 0];
    etaTarget=[etaTarget 1];
    
    outname = sprintf('./data/D11%i_datat.data',macro_meso_iteration);
    csvwrite(outname,D11);
     
    outname = sprintf('./data/D12%i_datat.data',macro_meso_iteration);
     csvwrite(outname,D12);
     
       outname = sprintf('./data/D22%i_datat.data',macro_meso_iteration);
     csvwrite(outname,D22);
     
        outname = sprintf('./data/D33%i_datat.data',macro_meso_iteration);
     csvwrite(outname,D33);
     
     % write the targets
     outname = sprintf('./data/pstrain1%i_datat.data',macro_meso_iteration);
     csvwrite(outname,pstrain1);
     
    outname = sprintf('./data/pstrain2%i_datat.data',macro_meso_iteration);
     csvwrite(outname,pstrain2);
     
       outname = sprintf('./data/pstrain3%i_datat.data',macro_meso_iteration);
     csvwrite(outname,pstrain3);
     
        outname = sprintf('./data/etaTarget%i_datat.data',macro_meso_iteration);
        csvwrite(outname,etaTarget);
        
elseif(step==3)
      macro_meso_iteration=1;
    outname = sprintf('./data/D11%i_datat.data',macro_meso_iteration);
    D11=  csvread(outname);
    
    outname = sprintf('./data/D12%i_datat.data',macro_meso_iteration);
    D12=  csvread(outname);
    
    outname = sprintf('./data/D22%i_datat.data',macro_meso_iteration);
    D22=  csvread(outname);
    
    outname = sprintf('./data/D33%i_datat.data',macro_meso_iteration);
    D33=  csvread(outname);
    
    % write the targets
    outname = sprintf('./data/pstrain1%i_datat.data',macro_meso_iteration);
    pstrain1=  csvread(outname);
    
    outname = sprintf('./data/pstrain2%i_datat.data',macro_meso_iteration);
    pstrain2=  csvread(outname);
    
    outname = sprintf('./data/pstrain3%i_datat.data',macro_meso_iteration);
    pstrain3=  csvread(outname);
    
    outname = sprintf('./data/etaTarget%i_datat.data',macro_meso_iteration);
    etaTarget=  csvread(outname);
    
    
    %-----------------------------
    % Plot data
    %-----------------------------
    plotting = 1;
    if plotting ==1
        maximizePlots=1;
        
        %     RhoColumn=D12; % color
        %     circleSize = ones(size(RhoColumn))*100; % circle size.
        %     scatter3(D11,D22,D33,circleSize, pstrain1,'filled','MarkerEdgeColor','k')
        %     title(sprintf('d11 plot as a function of e1 e2, e3'));
        %     colorbar
        %     xlabel('D11');
        %     ylabel('D22');
        %     zlabel('D33');
        %     colormap winter
        D_array=[D11;D12;D22;D33];
        
        names = {'d11' 'd12' 'd22' 'd33'}
        count=1;
        for i = 1:4
            count=1;
            for dd = densityTargetsVector
                subplot(5,4,count)
                
                targetDensity=dd ;
                count=count+1;
                error=0.01;
                logicarray = etaTarget>targetDensity-error;
                logicarray2 = etaTarget<targetDensity+error;
                logic3 = logicarray+logicarray2;
                
                pstrain1_temp=pstrain1(logic3==2);
                pstrain2_temp=pstrain2(logic3==2);
                pstrain3_temp=pstrain3(logic3==2);
                
                Dvalues = D_array(i,:);
                D_indexPosition=Dvalues(logic3==2);
                
                RhoColumn=D_indexPosition; % color
                circleSize = ones(size(RhoColumn))*100; % circle size.
                
                
                scatter3(pstrain1_temp,pstrain2_temp,pstrain3_temp,circleSize, RhoColumn,'filled','MarkerEdgeColor','k')
                title(sprintf('%s plot: density = %f',char(names(i)),targetDensity));
                colorbar
                xlabel('pstrain1');
                ylabel('pstrain2');
                zlabel('pstrain3');
                colormap winter
%                 caxis([0 matProp.E_material1 ]);
                %colormap('gray')
                % colormap(flipud(gray(256)));
                % colormap('summer');
                
                
            end
            if(maximizePlots==1)
                set(gcf, 'Position', get(0, 'Screensize'));
            end
            nameGraph2 = sprintf('./data/%sfunctionOfPseudoStrains%i.png', char(names(i)),config.macro_meso_iteration);
            print(nameGraph2,'-dpng');
            close all
        end
    else
        
        
%         scatter3(D11,D22,D33,circleSize, pstrain1)
%         logicarray = pstrain1>-10;
%         logicarray2 =  pstrain2>0;
%         logic3 = logicarray+1;
%         pstrain1_temp=pstrain1(logic3==2);
%         pstrain2_temp=pstrain2(logic3==2);
%         pstrain3_temp=pstrain3(logic3==2);
%         
%         etaTarget_temp=etaTarget(logic3==2);
        
        D11_temp=D11;
        D12_temp=D12;%(logic3==2);
        D22_temp=D22;%(logic3==2);
        D33_temp=D33;%(logic3==2);
        
        
        % x = [D11 ;  D22; D33 ];
        x = [D11_temp ;  D22_temp; D33_temp ];
        %     t=[pstrain1;pstrain2;pstrain3; etaTarget ];
        % t=[etaTarget_temp;pstrain1_temp;pstrain2_temp;pstrain3_temp];
        % t=[etaTarget_temp;pstrain1_temp;pstrain2_temp;pstrain3_temp];
        t=[pstrain1_temp];
        size(x)
        size(t)
        
        %         setdemorandstream(491218382)
        % net = fitnet(10,'trainbfg');
        %           net = cascadeforwardnet(100);
        net = fitnet(40,'trainbfg');
        %         net = feedforwardnet(80);
        figure(1)
        % view(net)
        % print('ANN_view.png','-dpng');
        
        % Train
        [net,tr] = train(net,x,t);
        tr
        
        plotperform(tr)
        print('ANN_preformance.png','-dpng');
        
        
        y = net(x);
        plotregression(t,y)
        print('ANN_regressionTest2.png','-dpng');
        
        genFunction(net,'predictPseudoStrainAndTarget.m');
    end
end





    


%ExxValues=padarray(ExxValues,
% ExxValues=reshape(ExxValues,1,[]);

end