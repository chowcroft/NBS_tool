function O = rigidPart_ExampleHub(master_object,parent_object,name)

O = NBS_rigidPart(master_object,name,'Parent',parent_object);

%--------------------------------------------------------------------------
L = 2;
O.s = reshape(linspace(-L,L/2,50),1,1,[]); ns = numel(O.s);
O.s_aero = reshape([0 L],1,1,[]);
O.isAero = false;
    
O.Gamma_W = O.s.*[0;1;0];

temp = linspace(0,1,ns);
O.w = reshape(temp.^0.3,1,1,[])*2;
O.h = O.w;
%--------------------------------------------------------------------------
O.root_idx = 1;
O.Mu = 10;
O.rMu_W = [0;L/2*0+1;0];
%--------------------------------------------------------------------------
I_tau = 1;
O.I_varTheta = [
    I_tau 0     0;
    0 I_tau 0;
    0 0     I_tau];
%--------------------------------------------------------------------------
O.CrossSectionProfiles = 'circ';
O.root_idx = 1;

end