function O = testCase_ClampedPatilHodgesWing_BeamPropertiesObject()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
%O.aerodynamics = 'strip_unsteady';
O.uVec_freeStream_G = [1;0;0];
%--------------------------------------------------------------------------
O.grav_acc = 9.807*0;
O.gravVec_G = [0;0;-1];

%--------------------------------------------------------------------------
O.qRigidT = [];
O.qRigidR = [];

%==========================================================================
%>>Add Flexible Patil/Hodges Half-Wing

O_halfWing = parameter_sets.FlexPart_Library.flexPart_Patil_Hodges_16m_halfWing_BeamPropertiesObject(O,O,'halfWing');
O_halfWing.wingRoot_offset_A = [0;0;0];

%//////Applied Loads
O_halfWing.tip_force_global  = [100;0;100];
O_halfWing.tip_force_local   = [0;0;0];
O_halfWing.tip_moment_global = [0;0;0];
O_halfWing.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_halfWing.shape_class_bend = 'chebyshev_1st';
O_halfWing.shape_class_twist = 'chebyshev_1st';
O_halfWing.shape_class_shear = 'chebyshev_1st';
O_halfWing.shape_class_exten = 'chebyshev_1st';
O_halfWing.shape_BCs_bend  = [0 1;1 1;1 1];
O_halfWing.shape_BCs_twist = [0 1;1 1;1 1];
O_halfWing.shape_BCs_shear = [1 1;1 1;1 1];
O_halfWing.shape_BCs_exten = [1 1;1 1;1 1];

O_halfWing.qth.n = 5; O_halfWing.qth.group = 'qth1';
O_halfWing.qsi.n = 5; O_halfWing.qsi.group = 'qsi1';
O_halfWing.qph.n = 5; O_halfWing.qph.group = 'qph1';
O_halfWing.qSx.n = 0; O_halfWing.qSx.group = 'qSx1';
O_halfWing.qSy.n = 0; O_halfWing.qSy.group = 'qSy1';
O_halfWing.qSz.n = 0; O_halfWing.qSz.group = 'qSz1';

O_halfWing.populate_shape_set('PLOT',false);

% %------modify shape sets to allow for kink discontinuity-------------------
% sob = O_halfWing.shapeObject_bend;
% sot = O_halfWing.shapeObject_twist;
% 
% sob.addCustomFunctions(obj,s_custom,y_custom,dy_ds_custom,varargin)

%--------------------------------------------------------------------------

% O.flexParts_nonlinear.halfWing = O_halfWing;
% O.flexParts_nonlinear_cell = [O.flexParts_nonlinear_cell ; {O_halfWing}];
%==========================================================================

O.plotBounds = [[-0.5 0.5];[0 1];[-0.5 0.5]]*O_halfWing.L*1.1;
set_dependent_properties(O);

end