
coelhoMethod = 1;
DiffStore=[];
if(coelhoMethod==1)
    folderNum=0;
    for i = 1:20
        % get the density field
        outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,i);
        xnew = csvread(outname);
        
        if(i>1)
            DiffX = xold-xnew;
            sumOfDiffX = sum(sum(abs(DiffX)));
            
            %if(sumOfDiffX<=0.005)
            fprintf('%i,Diff = %f\n',i,sumOfDiffX)
            DiffStore=[DiffStore sumOfDiffX];
            %end
        end
        
        xold = xnew;
        
    end
    plot(DiffStore)
    title('Summed Difference in density values')
    xlabel('Iterations')
       ylabel('Summed Diff')
      %  nameGraph = sprintf('./gradTopOptimization%fwithmesh%i_load%i.png', config.w1,config.macro_meso_iteration,loadcaseIndex);
            print('SummedDiffPlotCoelhoMethod','-dpng');
end

