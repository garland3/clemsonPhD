function   [ checkedElements] = CalculateCheckedElements(ne, settings)

count = 1;
checkedElements=[];

evenOdd = mod(settings.macro_meso_iteration,2);

for i = 1:settings.nelx
    for j =1: settings.nely
         temp = mod(i +j,2);
      
         if(temp==evenOdd ||settings.macro_meso_iteration>4 )
             checkedElements=[checkedElements count];
         end
        count=count+1;
    end
end
