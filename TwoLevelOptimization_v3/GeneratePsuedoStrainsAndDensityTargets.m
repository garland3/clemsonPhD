function [] =   GeneratePsuedoStrainsAndDensityTargets(DV,config,matProp, step)

% Usize = size(config.loadingCase,2);
% nnodes=config.targetTestVectorLen^4*2;
% DV.U=ones(Usize,nnodes);% 1; % so it is not NULL. Just assign a placeholder value
% DV.sensitivityElastic=1;
% DV.sensitivityElasticPart2=1;
% DV.w=1;
% DV.lambda1=1;
method = 3;
totalNotInUse=1;
if (method ==1)
    temp=config.targetTestVectorLen-1;
    term1 =0:0.5/temp:0.5; % -1 to 1 is the domain
    term2 =0:0.5/temp:0.5;  % -1 to 1 is the domain
    term3 =-0.3:1/temp:0.3;  % to 1 is the domain
    densityTargetsVector = 2/temp:1/temp:1-1/temp;  % 0 to 1 is the domain
    [strain1, strain2, strain3, densityTargets] = ndgrid(term1,term2,term3,densityTargetsVector);
    
elseif(method==2)
    temp=config.targetTestVectorLen-1;
    term1 =0:1/temp:0.5; % -1 to 1 is the domain
    term2 =0:1/temp:0.5;  % -1 to 1 is the domain
    term3 =0:1/temp:0.3;  % to 1 is the domain
    densityTargetsVector = 2/temp:1/temp:1-1/temp;  % 0 to 1 is the domain
    [strain1_temp, strain2_temp, strain3_temp, densityTargets_temp] = ndgrid(term1,term2,term3,densityTargetsVector);
    
    % force B1 > B2 always
    logicArray = strain1_temp<strain2_temp;
    totalNotInUse = sum(sum(sum(sum(logicArray))));
    dontUseArray = ones(size(strain1_temp))*-1;
    %     strain1=strain1_temp(logicArray);
    %     strain2=strain2_temp(logicArray);
    %     strain3=strain3_temp(logicArray);
    %     densityTargets=densityTargets_temp(logicArray);
    strain1=strain1_temp;
    strain2=strain2_temp;
    strain3=strain3_temp;
    densityTargets=densityTargets_temp;
    strain1(logicArray)=dontUseArray(logicArray);
    strain2(logicArray)=dontUseArray(logicArray);
    strain3(logicArray)=dontUseArray(logicArray);
    densityTargets(logicArray)=dontUseArray(logicArray);
    
    
    % force ratio of B1/B2 or B2/B1 to not be larger than ratioLimit
    %     ratioLimit=2;
    %     logicArray1 =( ( strain2./strain1<ratioLimit)+ densityTargets>0.3);
    %     logicArray=logicArray1(logicArray1>=1);
    %     strain1=strain1(logicArray);
    %     strain2=strain2(logicArray);
    %     strain3=strain3(logicArray);
    %     densityTargets=densityTargets(logicArray);
    %
    
elseif(method==3)
    temp=config.targetTestVectorLen-1;
    term1 =0:1/temp:0.5; % -1 to 1 is the domain
    term2 =0:1/temp:0.5;  % -1 to 1 is the domain
    term3 =0:1/temp:0.3;  % to 1 is the domain
    densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
    [strain1_temp, strain2_temp, strain3_temp, densityTargets_temp] = ndgrid(term1,term2,term3,densityTargetsVector);
    
    % force B1 > B2 always, aCtually, we can check for this latter. We
    % don't want to mess up the ndgrid data.
    logicArray = strain1_temp<strain2_temp;
    totalNotInUse = sum(sum(sum(sum(logicArray))));
    %     dontUseArray = ones(size(strain1_temp))*-1;
    %     strain1=strain1_temp(logicArray);
    %     strain2=strain2_temp(logicArray);
    %     strain3=strain3_temp(logicArray);
    %     densityTargets=densityTargets_temp(logicArray);
    strain1=strain1_temp;
    strain2=strain2_temp;
    strain3=strain3_temp;
    densityTargets=densityTargets_temp;
    %     strain1(logicArray)=dontUseArray(logicArray);
    %     strain2(logicArray)=dontUseArray(logicArray);
    %     strain3(logicArray)=dontUseArray(logicArray);
    %     densityTargets(logicArray)=dontUseArray(logicArray);
    
end

SIMPTemp = reshape(strain1,[],1);
[ numProblems t1]=size(SIMPTemp);
numProblemsActual=numProblems-totalNotInUse

% strain1_before = strain1;
strain1 = reshape(strain1,[],1);
strain2 = reshape(strain2,[],1);
strain3 = reshape(strain3,[],1);


t1=size(term1,2);
t2=size(term2,2);
t3=size(term3,2);
t4=size(densityTargetsVector,2);
% strain1_after = reshape(strain1,[t1 t2 t3 t4]);

%  logicArray = abs(strain1_before-strain1_after)<0.00001




densityTargets = reshape(densityTargets,[],1);
folderNum=0;
if(step==1)
    fprintf('Generate target values psuedo strain and density target test\n');
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
    
    D11Local = 0;
    D12Local= 0;
    D22Local = 0;
    D33Local = 0;
    
    for e = 1: ne
        fprintf('element %i of %i\n',e,ne);
        
        p1_local = strain1FromCSV(e);
        p2_local =strain2FromCSV(e);
        p3_local = strain3FromCSV(e);
        etaLocal = densityTargetFromCSV(e);
        
        if(p1_local>=p2_local)
            %             fprintf('Do not use this set\n')
            %             continue
            %         end
            
            
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
            
            
            %   fprintf('D11 \t\t%f, \nD22 \t\t%f\nD33 \t\t%f\npstrain \t%f %f %f\nDensity \t%f\n\n',D11Local,D11Local,D33Local,P1_local,p2_local,p3_local,etaLocal);
            
            outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,e);
            etaArray = csvread(outname);
            
            if 1==0
                p.PlotArrayGeneric(etaArray,'eta')
                title({sprintf('pstrain %f %f %f',p1_local,p2_local,p3_local),sprintf('Density %f',etaLocal),sprintf('D11 %f, D22 %f,D33 %f',D11Local,D11Local,D33Local )  })
                caxis([0 1])
                drawnow
                
            end
        else
            % Just use the previous last known value.
            %              D11Local = 0;
            %             D12Local= 0;
            %             D22Local = 0;
            %             D33Local = 0;
            
        end
        
        D11=[D11; D11Local];
        D12=[D12; D12Local];
        D22=[D22; D22Local];
        D33=[D33; D33Local];
        
        pstrain1=[pstrain1; p1_local];
        pstrain2=[pstrain2; p2_local];
        pstrain3=[pstrain3; p3_local];
        etaTarget=[etaTarget; etaLocal];
        
        
    end
    
    %     % ---------------
    %     % Add a zero case
    %     %-----------------
    %     D11=[D11 0];
    %     D12=[D12 0];
    %     D22=[D22 0];
    %     D33=[D33 0];
    %
    %     pstrain1=[pstrain1 0];
    %     pstrain2=[pstrain2 0];
    %     pstrain3=[pstrain3 0];
    %     etaTarget=[etaTarget 0];
    %
    %      % ---------------
    %     % Add a one case case , totally solid
    %     %-----------------
    %     D11=[D11 matProp.E_material1/(1-matProp.v^2)];
    %     D12=[D12 matProp.v*matProp.E_material1/(1-matProp.v^2)];
    %     D22=[D22 matProp.E_material1/(1-matProp.v^2)];
    %     D33=[D33 1/2*(1-matProp.v)*matProp.E_material1/(1-matProp.v^2)];
    %
    %     pstrain1=[pstrain1 0];
    %     pstrain2=[pstrain2 0];
    %     pstrain3=[pstrain3 0];
    %     etaTarget=[etaTarget 1];
    
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
    
    % read the targets
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
    runANN=0;
    showAllDesigns =1;
    if plotting ==1
        if(method==2 || method==3)
            % Remove the values that don't meet the logic test.
            logticTest=pstrain1>=pstrain2;
            D11=D11(logticTest);
            D12=D12(logticTest);
            D22=D22(logticTest);
            D33=D33(logticTest);
            pstrain1=pstrain1(logticTest);
            pstrain2=pstrain2(logticTest);
            pstrain3=pstrain3(logticTest);
            etaTarget=etaTarget(logticTest);
            
        end
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
        
        makevide=1;
        if(makevide ==1)
            config.recvid=1;
            video = VideoManager;
            [vidObj, framedNumber] = video.InitializeVideo( config,'./psuedoStrainDensityTarget.avi');
            F=getframe();
        end
        
        count=1;
        D_array=[D11 D12 D22 D33];
        
        names = {'d11' 'd12' 'd22' 'd33'}
        count=1;
        llArray =[1 2 5 9];
        for dd = densityTargetsVector
            for i = 1:4
                
                ll = llArray(i);
                %                                  subplot(3,3,count)
                subplot(3,3,ll)
                
                targetDensity=dd ;
                count=count+1;
                error=0.01; % needed becasue of numeric issues and 1~=1 when rounded off.
                logicarray = etaTarget>targetDensity-error;
                logicarray2 = etaTarget<targetDensity+error;
                logic3 = logicarray+logicarray2;
                
                pstrain1_temp=pstrain1(logic3==2);
                pstrain2_temp=pstrain2(logic3==2);
                pstrain3_temp=pstrain3(logic3==2);
                
                Dvalues = D_array(:,i);
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
                
                
                %end
                if(maximizePlots==1)
                    set(gcf, 'Position', get(0, 'Screensize'));
                end
                %             nameGraph2 = sprintf('./data/%sfunctionOfPseudoStrains%i.png', char(names(i)),config.macro_meso_iteration);
                
                
            end
            nameGraph2 = sprintf('./data/plots/%sfunctionOfPseudoStrains%i_%i.png', char(names(i)),config.macro_meso_iteration,count);
            print(nameGraph2,'-dpng');
            count = count+1;
            
            if(makevide ==1)
                [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
            end
            
            close all
        end
        
        if(makevide ==1)
            video.CloseVideo( config, F,vidObj)
        end
        %
    end
    if(runANN==1)
        if(method==2 || method==3)
            % Remove the values that don't meet the logic test.
            logticTest=pstrain1>=pstrain2;
            D11=D11(logticTest);
            D12=D12(logticTest);
            D22=D22(logticTest);
            D33=D33(logticTest);
            pstrain1=pstrain1(logticTest);
            pstrain2=pstrain2(logticTest);
            pstrain3=pstrain3(logticTest);
            etaTarget=etaTarget(logticTest);
            
        end
        
        
        %         scatter3(D11,D22,D33,circleSize, pstrain1)
        logicarray = pstrain1>=pstrain2;
        % logicarray2 =  pstrain2>0;
        % logic3 = logicarray+1;
        pstrain1_temp=pstrain1(logicarray==1);
        pstrain2_temp=pstrain2(logicarray==1);
        pstrain3_temp=pstrain3(logicarray==1);
        
        etaTarget_temp=etaTarget(logicarray==1);
        
        D11_temp=D11(logicarray==1);
        D12_temp=D12(logicarray==1);%(logic3==2);
        D22_temp=D22(logicarray==1);%(logic3==2);
        D33_temp=D33(logicarray==1);%(logic3==2);
        
        
        % x = [D11 ;  D22; D33 ];
        x = [D11_temp  D12_temp D22_temp D33_temp ];
        %     t=[pstrain1;pstrain2;pstrain3; etaTarget ];
        % t=[etaTarget_temp;pstrain1_temp;pstrain2_temp;pstrain3_temp];
        t=[etaTarget_temp pstrain1_temp pstrain2_temp pstrain3_temp];
        %         t=[pstrain1_temp ; pstrain2_temp];
        size(x)
        size(t)
        
        useParallel = 1;
        if(useParallel==1)
            poolobj = gcp('nocreate'); % If no pool,create new one.
            if isempty(poolobj)
                parpool('local', 3)
                poolsize = 3;
            else
                poolsize = poolobj.NumWorkers;
            end
            poolsize
        end
        
        setdemorandstream(491218382)
        %          net = fitnet(10,'trainbfg');
        %  net = cascadeforwardnet(10);
        %         net = feedforwardnet(3);
        net = feedforwardnet(10,'trainscg');
        figure(1)
        % view(net)
        % print('ANN_view.png','-dpng');
        
        % Train
        if(useParallel~=1)
            [net,tr] = train(net,x,t);
            y = net(x);
        else
            [net,tr] = train(net,x,t,'useParallel','yes')
            y = net(x,'useParallel','yes');
        end
        
        plotperform(tr)
        print('ANN_preformance.png','-dpng');
        
        plotregression(t,y)
        print('ANN_regressionTest2.png','-dpng');
    end
    if(showAllDesigns==1)
        % -----------------------------
        % -- PLOT the designs
        % --------------------------
        
        %         temp=config.targetTestVectorLen-1;
        %         term1 =0:1/temp:0.5; % -1 to 1 is the domain
        %         term2 =0:1/temp:0.5;  % -1 to 1 is the domain
        %         term3 =0:1/temp:0.3;  % to 1 is the domain
        %         densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
        %
        %         t1=size(term1,2);
        %         t2=size(term2,2);
        %         t3=size(term3,2);
        %         t4=size(densityTargetsVector,2);
        
        %                 D11=D11(1:end-2);
        %                 D12=D12(1:end-2);
        %                 D22=D22(1:end-2);
        %                 D33=D33(1:end-2);
        %                 etaTarget=etaTarget(1:end-2);
        %                 pstrain1=pstrain1(1:end-2);
        %                 pstrain2=pstrain2(1:end-2);
        %                 pstrain3=pstrain3(1:end-2);
        elementNumbers = 1:t1*t2*t3*t4;
        
        
        %         D11 = reshape(D11,[t1 t2 t3 t4]);
        %         D12 = reshape(D12,[t1 t2 t3 t4]);
        %         D22 = reshape(D22,[t1 t2 t3 t4]);
        %         D33 = reshape(D33,[t1 t2 t3 t4]);
        elementNumbers=reshape(elementNumbers,[t1 t2 t3 t4]);
        
        
        
        etaTarget = reshape(etaTarget,[t1 t2 t3 t4]);
        pstrain1 = reshape(pstrain1,[t1 t2 t3 t4]);
        pstrain2 = reshape(pstrain2,[t1 t2 t3 t4]);
        pstrain3 = reshape(pstrain3,[t1 t2 t3 t4]);
        
        
        makevide=1;
        if(makevide ==1)
            config.recvid=1;
            video = VideoManager;
            [vidObj, framedNumber] = video.InitializeVideo( config,'./valid/mesoValidation.avi');
            F=getframe();
        end
        
        config.   numTilesX = 1;
        config.     numTilesY = 1;
        
        config.nelx=t1;
        config.nely=t2;
        
        configMeso = config;
        configMeso.nelx=configMeso.nelxMeso;
        configMeso.nely=configMeso.nelxMeso;
        
        term1 =0:1/temp:0.5; % -1 to 1 is the domain
        term2 =0:1/temp:0.5;  % -1 to 1 is the domain
        term3 =0:1/temp:0.3;  % to 1 is the domain
        densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
        %     [strain1_temp, strain2_temp, strain3_temp, densityTargets_temp] = ndgrid(term1,term2,term3,densityTargetsVector);
        
        t1=size(term1,2);
        t2=size(term2,2);
        t3=size(term3,2);
        t4=size(densityTargetsVector,2);
        
        p = plotResults;
        % plot each eta on a different plot, and make a video.
        %         elementCount=1;
        figure
        for i = 1:t4
            % plot each p3 on a different plot
            for j = 1:t3
                % make a blank huge array
                totalY = t1*configMeso.nelxMeso;
                completeStruct = zeros(totalY,totalY);
                xxx=ones(t1,t2);
                for k = 1:t2
                    for l = 1:t1
                        eta = etaTarget(l,k,j,i);
                      
                        
                        elementCount=elementNumbers(l,k,j,i);
                        p3=pstrain3(l,k,j,i);
                        fprintf('element %i  with eta %f and ps(3) %f\n',elementCount,eta,p3);
                        
                        if(  pstrain2(l,k,j,i)> pstrain1(l,k,j,i) )
                            %                             macroEleProps.targetDensity=-1;
                            fprintf('Do not use in test, ps2>ps1 symmetric\n');
                            continue;
                        end
                        
                        macroElementProps.elementNumber=elementCount;
                        %                 results = elementXYposition(macroElementProps.elementNumber,:);
                        macroElementProps.yPos = k;
                        macroElementProps.xPos = l;
                        
                        x=GetMesoUnitCellDesignFromCSV(config,elementCount);
                        DV.x = x;
                        
                       
                        
                        step = 1;
                        completeStruct= TileMesoStructureV2(configMeso,config, DV,macroElementProps,xxx,completeStruct,step);
                        %                         elementCount=elementCount+1;
                    end
                end
                completeStruct( completeStruct>1)=1;
                
                completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
                completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;                
        
                plotname = sprintf('meso designs  with eta %f and ps(3) %f',eta,p3);
                p.PlotArrayGeneric( completeStruct, plotname)
                
                xlabel(gca,'p1');
                ylabel(gca,'p2');
                set(get(gca,'YLabel'),'visible','on') % not sure why this is needed. 
                set(get(gca,'XLabel'),'visible','on')
                
                drawnow               
             
                nameGraph = sprintf('./valid/validationPlot_etaIndex%i_psIndex%i.png', i,j);
                %         print(nameGraph,'-dpng', '-r1200')
                print(nameGraph,'-dpng')                
                
                if(makevide ==1)
                    [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
                end
                
                
            end
        end
        if(makevide ==1)
            video.CloseVideo( config, F,vidObj)
        end
    end
    
    
    
    
    
    
    
    
    %ExxValues=padarray(ExxValues,
    % ExxValues=reshape(ExxValues,1,[]);
    
end