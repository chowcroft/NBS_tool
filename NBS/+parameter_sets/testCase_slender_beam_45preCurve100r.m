function O = testCase_slender_beam_45preCurve100r()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
O.aerodynamics = '';
O.uVec_freeStream_G = [1;0;0];
%--------------------------------------------------------------------------
O.grav_acc = 9.807*0;
O.gravVec_G = [0;0;-1];
%--------------------------------------------------------------------------
O.qRigidT = [];
O.qRigidR = [];
%==========================================================================

%==========================================================================
%>>Add Part

O_slenderBeamPreCurve = parameter_sets.FlexPart_Library.flexPart_slender_beam_45preCurve100r(O,O,'slender_beam_45preCurve100r');
O_slenderBeamPreCurve.wingRoot_offset_A = [0;0;0];

%//////Applied Loads
O_slenderBeamPreCurve.tip_force_global  = [0;0;0];
O_slenderBeamPreCurve.tip_force_local   = [0;0;600];
O_slenderBeamPreCurve.tip_moment_global = [0;0;0];
O_slenderBeamPreCurve.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_slenderBeamPreCurve.shape_class_bend = 'chebyshev_1st';
O_slenderBeamPreCurve.shape_class_twist = 'chebyshev_1st';
O_slenderBeamPreCurve.shape_class_shear = 'chebyshev_1st';
O_slenderBeamPreCurve.shape_class_exten = 'chebyshev_1st';
O_slenderBeamPreCurve.shape_BCs_bend  = [0 1;1 1;1 1];
O_slenderBeamPreCurve.shape_BCs_twist = [0 1;1 1;1 1];
O_slenderBeamPreCurve.shape_BCs_shear = [1 1;1 1;1 1];
O_slenderBeamPreCurve.shape_BCs_exten = [1 1;1 1;1 1];
%//////State Assignment
O_slenderBeamPreCurve.qth.n = 8; O_slenderBeamPreCurve.qth.group = 'qth1';
O_slenderBeamPreCurve.qsi.n = 8; O_slenderBeamPreCurve.qsi.group = 'qsi1';
O_slenderBeamPreCurve.qph.n = 8; O_slenderBeamPreCurve.qph.group = 'qph1';
O_slenderBeamPreCurve.qSx.n = 0; O_slenderBeamPreCurve.qSx.group = 'qSx1';
O_slenderBeamPreCurve.qSy.n = 0; O_slenderBeamPreCurve.qSy.group = 'qSy1';
O_slenderBeamPreCurve.qSz.n = 0; O_slenderBeamPreCurve.qSz.group = 'qSz1';

O_slenderBeamPreCurve.populate_shape_set('PLOT',false);
%==========================================================================

O.plotBounds = [-30 30;-40 80;-10 70];
set_dependent_properties(O);

end