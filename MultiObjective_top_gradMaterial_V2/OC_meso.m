%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc ,designVar, settings)

% dc = -dc; % maximization problem. 

l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
     lmid = 0.5*(l2+l1);
%     beta = -x.*sqrt(dc./lmid);
    beta = (sqrt(dc./lmid));
    lmid = 0.5*(l2+l1);
    xnew = max(0.01, ...
             max(x-move, ...
                 min(1., ....
                     min(x+move,x.*beta) ...
                    ) ...
                )...
              );
    
    %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
    %[volume1, volume2] = designVar.CalculateVolumeFractions(settings);
    %currentvolume=volume1+volume2;
    %if currentvolume - volfrac > 0;
    
    if sum(sum(xnew)) - volfrac*nelx*nely > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
end