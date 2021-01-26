function O = flexPart_Patil_Hodges_16m_halfWing(master_object,parent_object,name)

O = NBS_flexPart_nonlinear(master_object,name,'Parent',parent_object);

%--------------------------------------------------------------------------
L = 16;
ns = 201;
O.s = permute(linspace(0,L,ns),[1 3 2]);

nAnodes = 17;
Anode_Skew_Factor = 1.0; %dictates the degree to which the aero nodes are bunched towards the tip (set to 1 for linear distribution)
Anode_distr = linspace(1,0,nAnodes).^Anode_Skew_Factor;
s_aero_ = (1 - Anode_distr)*L;
O.s_aero = permute(s_aero_,[1 3 2]);
O.isAero = true;

O.w = ones(1,1,ns)*1;
O.h = ones(1,1,ns)*0.2;
%--------------------------------------------------------------------------
EIxx = 2e4; EIzz = 4e6; GJ  = 1e4;
O.StiffnessMatrix = [
    EIxx 0    0;
    0    GJ   0;
    0    0    EIzz];
O.StiffnessMatrix(4:6,4:6) = eye(3)*1e3;
%--------------------------------------------------------------------------
damping_factor = 0.04*100;
O.DampingMatrix = damping_factor*O.StiffnessMatrix;
%--------------------------------------------------------------------------
I_tau = 0.1;
O.I_varTheta_ps_I = [
    I_tau 0 0;
    0 I_tau 0;
    0 0 I_tau];
O.I_varTheta_discrete_I = 0;
%--------------------------------------------------------------------------
O.mps = 0.75;
O.msDiscrete = 0;
%--------------------------------------------------------------------------
O.R_W_WE = eye(3);
O.KAPPA_0_I = zeros(3,1,ns);
%--------------------------------------------------------------------------
O.c = ones(1,1,ns);
%O.aero_cntr = 0.25;
O.beam_cntr = ones(1,1,ns)*0.5;
%--------------------------------------------------------------------------
O.AICs = sqrt(1-(0.5*(O.s_aero(2:end)+O.s_aero(1:end-1))).^2./L^2)*1+0;
O.CrossSectionProfiles = 'NACA0012';
%--------------------------------------------------------------------------

end

