%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc ,designVar, settings)

% if the sensitivity is really small, then make it larger to help the
% optimal criteria method work better.
absL =  max(max(abs(dc)));
if absL <1000
    dc = dc*100000.0/absL;
end

multiplier = 1;

if(settings.doUseMultiElePerDV) % if elements per design var.     
   multiplier = settings.numVarsX*settings.numVarsY;   
else
    multiplier=nelx*nely;
end

l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    xnew = max(0.01,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
    
    %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
    
    %[volume1, volume2] = designVar.CalculateVolumeFractions(settings);
    %currentvolume=volume1+volume2;
    
    %if currentvolume - volfrac > 0;
    
    if sum(sum(xnew)) - volfrac*multiplier > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
end