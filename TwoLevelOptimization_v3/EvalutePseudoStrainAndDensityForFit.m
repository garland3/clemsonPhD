function diffOutput= EvalutePseudoStrainAndDensityForFit(x,D11sys,D12sys,D22sys,D33sys, ...
    D11_table4D,D12_table4D,D22_table4D,D33_table4D,...
    ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D)
% x(1) = ps(1) to evalute
% x(2) = ps(2) to evalute
% x(3) = ps(3) to evalute
% x(4) = etaTarget to evalute

% D11sys system value that we are trying to meet.
% D12sys system value that we are trying to meet.
% D22sys system value that we are trying to meet.
% D33sys system value that we are trying to meet.

% D11_table4D is the 4D array of known D11 values
% D12_table4D is the 4D array of known D12 values
% D22_table4D is the 4D array of known D22 values
% D33_table4D is the 4D array of known D33 values

% ps1_table4D is the 4D array of known ps values that coorospond with the D11_table4D
% ps2_table4D is the 4D array of known ps values that coorospond with the D11_table4D
% ps3_table4D is the 4D array of known ps values that coorospond with the D11_table4D
% eta_table4D is the 4D array of known eta values that coorospond with the D11_table4D

ps_new(1)=x(1);
ps_new(2)=x(2);
ps_new(3)=x(3);
eta_new=x(4);

if(ps_new(2)>  ps_new(1))
    temp = ps_new(2);
    ps_new(2)=ps_new(1);
    ps_new(1)=temp;
end


[D11_interp,D12_interp ,D22_interp,D33_interp]= InterpolatePseudoStrainsAndDensity(ps_new,eta_new, ...
    D11_table4D,D12_table4D,D22_table4D,D33_table4D,...
    ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D);



if(eta_new<-0.01)
    diffOutput=1e11;
    return
end

diffOutput = (D11_interp-D11sys)^2+(D12_interp-D12sys)^2+(D22_interp-D22sys)^2+(D33_interp-D33sys)^2;

