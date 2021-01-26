function O = flexPart_Patil_Hodges_16m_halfWing_BeamPropertiesObject(master_object,parent_object,partName)

%////////////////Structural Definition///////////////////

%create a blank beam properties object
BPObj = BeamPropertiesObject();

%populate the beam properties object with data at whatever resolution is supplied
L = 16;
ns = 201;
s_nodes = permute(linspace(0,L,ns),[1 3 2]);

BPObj.BeamTriad_reference = repmat(eye(3),1,1,ns);
BPObj.tangent_direction = 2;
BPObj.sBT = s_nodes;

EIxx = 2e4; EIzz = 4e6; GJ  = 1e4;
StiffMat = [
    EIxx 0    0;
    0    GJ   0;
    0    0    EIzz];
StiffMat(4:6,4:6) = eye(3)*1e3;
BPObj.StiffnessMatrix = repmat(StiffMat,1,1,ns);
BPObj.sK = s_nodes;
%BPObj.StiffnessFnc = @LinearStiffLaw <- default constitutive law

damping_factor = 0.2;
BPObj.DampingMatrix = damping_factor.*BPObj.StiffnessMatrix;
BPObj.sC = s_nodes;
%BPObj.DampingFnc = @LinearDampLaw <- default damping law

%BPObj.ReferenceOffset = default
%BPObj.sRO = default

I_tau = 0.1;
I_varTheta_ps_I = repmat(eye(3).*I_tau,1,1,ns);
mps = repmat(0.75,1,1,ns);

%BPObj.Mass_Continuous = default
%BPObj.Inertia_Continuous = default
%BPObj.MassOffset_Continuous = default
%BPObj.sM_C = default

BPObj.MassPerSpan_Continuous = mps;
BPObj.InertiaPerSpan_Continuous = I_varTheta_ps_I;
BPObj.MassPerSpanOffset_Continuous = zeros(3,1,ns);
BPObj.sMPS_C = s_nodes;

%BPObj.Mass_Discrete = default
%BPObj.Inertia_Discrete = default
%BPObj.MassOffset_Discrete = default
%BPObj.sM_D = default


%////////////////Aero Definition///////////////////
nAnodes = 21;
Anode_Skew_Factor = 1.0; %dictates the degree to which the aero nodes are bunched towards the tip (set to 1 for linear distribution)
Anode_distr = linspace(1,0,nAnodes).^Anode_Skew_Factor;
s_aero_ = (1 - Anode_distr)*L;
s_aero = permute(s_aero_,[1 3 2]);
%BPObj.alpha_A = decide on parameter definition for calculation of EAp_A
%BPObj.s_alpha = s_aero;

BPObj.chord = ones(1,1,nAnodes);
BPObj.s_chord = s_aero;

BPObj.beam_cntr = 0.5*ones(1,1,nAnodes);
BPObj.s_bcntr = s_aero;


%////////////////Create an NBS_flexPart_nonlinear object///////////////////

%combine the mass sources into a single discrete distribution
BPObj.unify_mass_distributions();
BPObj.aggregate_mass_sources();

%TODO make sure code can execute if mass sample distribution is different to 's'

%resample structural and aero parameters over common sample distributions
BPObj.change_samplePoints('structural',s_nodes);
BPObj.change_samplePoints('aero',s_nodes);

O = BPObj.NBS_export(master_object,partName,'Parent',parent_object);

%--------------------------------------------------------------------------
O.sweep_root = 0;%deg
O.dihedral_root = 0;%deg
O.alpha_root = 5;%deg
%update_R_A_W(O,'sweep_dihedral_alpha_deg',{sweep_root,dihedral_root,alpha_root},'reflect',false);
%O.E0_A = MultiProd_(O.R_A_W,O.E0_A);
%--------------------------------------------------------------------------
O.w = ones(1,1,ns)*1;
O.h = ones(1,1,ns)*0.2;
%--------------------------------------------------------------------------
O.s_aero = s_aero;
O.isAero = true;
O.AICs = sqrt(1-(0.5*(O.s_aero(2:end)+O.s_aero(1:end-1))).^2./L^2)*0+1;
O.CrossSectionProfiles = 'NACA0012';
%--------------------------------------------------------------------------

end

