function [D11_interp,D12_interp ,D22_interp,D33_interp]= InterpolatePseudoStrainsAndDensity(ps_new,eta_new, ...
    D11_table4D,D12_table4D,D22_table4D,D33_table4D,...
    ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D)


D11_interp = interpn(ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D,D11_table4D,ps_new(1),ps_new(2),ps_new(3),eta_new,'linear');
D12_interp = interpn(ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D,D12_table4D,ps_new(1),ps_new(2),ps_new(3),eta_new,'linear');
D22_interp = interpn(ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D,D22_table4D,ps_new(1),ps_new(2),ps_new(3),eta_new,'linear');
D33_interp = interpn(ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D,D33_table4D,ps_new(1),ps_new(2),ps_new(3),eta_new,'linear');
