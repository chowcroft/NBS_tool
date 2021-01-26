function O = testCase_temp_nblade_example()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master('partName','Example_Turbine');
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30/2; %m/s                                                             free stream airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
O.aerodynamics = 'strip_steady'; %                                         aerodynamic model
O.uVec_freeStream_G = [1;0;0]; %                                           free stream unit velocity vector
%--------------------------------------------------------------------------
O.grav_acc = 9.807*0;
O.gravVec_G = [0;0;-1];
%--------------------------------------------------------------------------
O.qRigidT = [0;0;0];
O.qRigidR = [0;0*pi/180;0];
%--------------------------------------------------------------------------
O.plotBounds = [[-1 1];[-1 1];[-1 1]]*40/2;

%==========================================================================

Rotor = struct();
nBlades = 4;

for i = 1:nBlades
    
    bladeName = ['Blade_', sprintf('%d',i)];
    O_blade = parameter_sets.FlexPart_Library.flexPart_Patil_Hodges_16m_halfWing(O,O,bladeName);
    
    radial_offset = 0.5;
    alpha_deg = 5;
    alpha_rad = alpha_deg*pi/180;
    
    method = 1;
    switch method
        case 1
            radialPosition_rad = 2*pi/nBlades*(i-1);
            O_blade.R_A_W = r_matrix([1,0,0],radialPosition_rad)*r_matrix([0,1,0],alpha_rad);
            O_blade.wingRoot_offset_A = O_blade.R_A_W*[0;radial_offset;0];
        case 2
            radialPosition_deg = 360/nBlades*(i-1)-1;
            O_blade.alpha_root = alpha_deg;
            O_blade.sweep_root = 0;
            if radialPosition_deg <= 90 || radialPosition_deg >= 270
                O_blade.dihedral_root = radialPosition_deg;
                O_blade.update_R_A_W();
            else
                O_blade.dihedral_root = 180-radialPosition_deg;
                O_blade.reflectedPart = true;
                O_blade.update_R_A_W();
            end
            O_blade.wingRoot_offset_A = O_blade.R_A_W*[0;radial_offset;0];
    end
    
    %//////Shape Functions
    O_blade.shape_class_bend = 'chebyshev_1st';
    O_blade.shape_class_twist = 'chebyshev_1st';
    O_blade.shape_BCs_bend  = [0 1;1 1;1 1];
    O_blade.shape_BCs_twist = [0 1;1 1;1 1];

    O_blade.qth.n = 8/2;
    O_blade.qsi.n = 6/2;
    O_blade.qph.n = 6/2;
    
    O_blade.populate_shape_set('PLOT',false,'setName',bladeName);
    
    currentBlade = ['O_turbineBlade_', num2str(i)];
    Rotor.(currentBlade) = O_blade;
    
end

%==========================================================================
%>>Add Hub

O_Hub = parameter_sets.RigidPart_Library.rigidPart_ExampleHub(O,O,'Hub');
O_Hub.connection_idx_ParentObj = 1;
O_Hub.connectionOffset_C = [0;0;0];
O_Hub.R_C_W_0 = r_matrix([0;0;-1],pi/2);

%==========================================================================
set_dependent_properties(O);
%==========================================================================

end