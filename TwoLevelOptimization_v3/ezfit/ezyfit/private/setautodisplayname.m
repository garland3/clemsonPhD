function setautodisplayname
%SETAUTODISPLAYNAME
%   SETAUTODISPLAYNAME gives an automatic 'DisplayName' field to the
%   data in the current plot which do not have already, in the form
%   'data1', 'data2'...
%
%   F. Moisy
%   Revision: 1.00,  Date: 2006/02/16
%
%   See also SHOWFIT, SELECTFIT.


% History:
% 2006/02/16: v1.00, first version.


ca=get(gca);
hcurve=ca.Children;
datanum=0;
for i=length(hcurve):-1:1,
    if isfield(get(hcurve(i)),'DisplayName')
        if isempty(get(hcurve(i),'DisplayName')),
            datanum = datanum+1;
            set(hcurve(i),'DisplayName',['data' num2str(datanum)]);
        end
    else
        return;
    end
end
