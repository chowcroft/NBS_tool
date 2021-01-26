function O = flexPart_slender_beam_45preCurve100r(master_object,parent_object,name)

O = NBS_flexPart_nonlinear(master_object,name,'Parent',parent_object);

%--------------------------------------------------------------------------
radius = 100;
L = 45*pi/180*radius;
ns = 101;
nAnodes = 21;

O.s = permute(linspace(0,L,ns),[1 3 2]);
O.s_aero = permute(linspace(0,L,nAnodes),[1 3 2]);
O.isAero = false;

O.w = 1;
O.h = 1;
density = 0.5;%4430/500;
%--------------------------------------------------------------------------
w = O.w; h = O.h;
A = w*h;
I1 = w*h^3/12;
I2 = w^3*h/12;
nn=[1:2:200]; xx8=tanh(nn*(pi*w/h/2))./(nn.^5); J=w*h^3/3*(1-192*h/w/pi^5*sum(xx8)); %given in pai book
E = 1e7;
poissonsRatio = 0*0.36;
G = E/(2*(1+poissonsRatio));
O.StiffnessMatrix(1:3,1:3) = [E*I1 0 0 ; 0 G*J 0 ; 0 0 E*I2];
O.StiffnessMatrix(4:6,4:6) = [A*G  0 0 ; 0 E*A 0 ; 0 0 A*G ];
%--------------------------------------------------------------------------
damping_factor = 1.0;
O.DampingMatrix = damping_factor*O.StiffnessMatrix;
%--------------------------------------------------------------------------
O.I_varTheta_ps_I = eye(3)*1e-3;
O.I_varTheta_ps_I(2,2) = density*(I1+I2);
O.I_varTheta_discrete_I = 0;
%--------------------------------------------------------------------------
O.mps = density*w*h;
O.msDiscrete = 0;
%--------------------------------------------------------------------------
%specification of curved shape
%=============
%difinition 1
O.KAPPA_0_I = -repmat([0;0;1]*(45*pi/180)/L,[1 1 ns]);
        O.th0 = +(45*pi/180)/L*O.s;
        O.dth_ds0 = +(45*pi/180)/L;
O.R_W_WE = [0 0 -1;0 1 0;1 0 0].';
%=============
%definition 2
% varTheta = reshape(linspace(0,-45,ns)*pi/180.*[0;0;1],3,1,[]);
% O.TD = rTransform_thetaVec_to_Rmat(varTheta);
% O.R_W_WE = eye(3);O.R_W_WE = [0 0 -1;0 1 0;1 0 0].';
%=============
%--------------------------------------------------------------------------
O.c = ones(1,1,ns);
% O.aero_cntr = 0.25; %<< currently the default
O.beam_cntr = ones(1,1,ns)*0.5;
%--------------------------------------------------------------------------
O.AICs = sqrt(1-(0.5*(O.s_aero(2:end)+O.s_aero(1:end-1))).^2./L^2);
O.CrossSectionProfiles = 'box';
%--------------------------------------------------------------------------
%O.nmo = O.Parent.nmo;
%--------------------------------------------------------------------------

end

