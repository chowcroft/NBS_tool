classdef BeamPropertiesObject < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
    
    properties%                                                            %[data structure][units] description
        
        BeamTriad_reference                                                %[(3)x(3)x(ns)][-] orthonormal intrinsic triad specifying the orientation of the beam
        tangent_direction                                                  %[{1,2,3}] specify which column of BeamTriad_reference runs tangent to the beam reference line
        sBT = []                                                           %[length ns vector][m] curvilinear sample points along the flexible beam
        %--                                                                %                                                                 |------curvature------| |-----shear-----|
        StiffnessMatrix                                                    %[(6)x(6)x(ns)][Nm^2] material stiffness matrix. Component order: kappa_x,kappa_y,kappa_z,tau_x,tau_y,tau_z. (x,y,z) coordinate system defined by the intrinsic triad 'obj.BeamTriad_reference'
        sK = []                                                            %[length ns vector][m] curvilinear sample points along the flexible beam
        StiffnessFnc = @LinearStiffLaw                                     %[function handle] Material stiffness function, defaults to linear relation
        %--
        DampingMatrix                                                      %[(6)x(6)x(ns)][Nm^2] material damping matrix. Component order: kappa_x,kappa_y,kappa_z,tau_x,tau_y,tau_z.
        sC = []                                                            %[length ns vector][m] curvilinear sample points along the flexible beam
        DampingFnc = @LinearDampLaw                                        %[function handle] Damping function, defaults to linear relation
        %--
        ReferenceOffset = 0                                                %[(3)x(2)x(ns)][m] reference line offsets for K and C matrices in the (x,y,z) 'obj.BeamTriad_reference' coordinate system. ReferenceOffset(:,1,:) <- left offsets , ReferenceOffset(:,2,:) <- right offsets.
        sRO = []                                                           %[length ns vector][m] curvilinear sample points along the flexible beam at which left/right offsets are given
        %--                                                                
        Mass_Continuous = 0                                                %[(1)x(1)x(ns)][kg] mass points (part of the continuous mass distribution)
        Inertia_Continuous = zeros(3)                                      %[(3)x(3)x(ns)][kg] inertia matrix, defined about 'obj.MassPerSpan' points and expressed in the 'obj.BeamTriad_reference' coordinate system
        MassOffset_Continuous = [0;0;0]                                    %[(3)x(1)x(ns)][m] mass offsets in the (x,y,z) 'obj.BeamTriad_reference' coordinate system
        sM_C = 0                                                           %[length ns vector][m] curvilinear sample points along the flexible beam
        %--                                                                
        MassPerSpan_Continuous = 0                                         %[(1)x(1)x(ns)][kg] mass points (part of the continuous mass distribution)
        InertiaPerSpan_Continuous = zeros(3)                               %[(3)x(3)x(ns)][kg] inertia matrix, defined about 'obj.MassPerSpan' points and expressed in the 'obj.BeamTriad_reference' coordinate system
        MassPerSpanOffset_Continuous = [0;0;0]                             %[(3)x(1)x(ns)][m] mass offsets in the (x,y,z) 'obj.BeamTriad_reference' coordinate system
        sMPS_C = 0                                                         %[length ns vector][m] curvilinear sample points along the flexible beam
        %--
        Mass_Discrete = 0                                                  %[(1)x(1)x(ns)][kg/m] discrete mass points
        Inertia_Discrete = zeros(3)                                        %[(3)x(3)x(ns)][kg/m] discrete inertia matrix, defined about 'obj.MassDiscrete' points and expressed in the 'obj.BeamTriad_reference' coordinate system
        MassOffset_Discrete = [0;0;0]                                      %[(3)x(1)x(ns)][m] mass offsets in the (x,y,z) 'obj.BeamTriad_reference' coordinate system
        sM_D = 0                                                           %[length ns vector][m] curvilinear sample points along the flexible beam
        %--
        
    end
    
    properties (Dependent)
        s_min
        s_max
        ReferenceRotation
        ReferenceTranslation
    end
    
    properties %Aerodynamic
        
        %alpha_A
        %s_alpha = []                                                      %[length ns vector][m] curvilinear sample points along the flexible beam
        %--
        chord
        s_chord                                                            %[length ns vector][m] curvilinear sample points along the flexible beam
        %--
        beam_cntr
        s_bcntr                                                            %[length ns vector][m] curvilinear sample points along the flexible beam
        %--
    end
    
    properties
        sampleVec_property_tags
    end
    
    methods %constructor method and NBS part export
        function obj = BeamPropertiesObject()
            obj.sampleVec_property_tags = {...
                'sBT'    ,'BeamTriad_reference'                            ,{'all','structural'};...
                'sK'     ,'StiffnessMatrix'                                ,{'all','structural','plot'};...
                'sC'     ,'DampingMatrix'                                  ,{'all','structural'};...
               %'s_alpha','alpha_A'                                        ,{'all','aero','plot'};...
                's_chord','chord'                                          ,{'all','aero','plot'};...
                's_bcntr','beam_cntr'                                      ,{'all','aero','plot'};...
                'sM_C'  ,'Mass_Continuous'                                 ,{'pointMass','plot'};...
                'sM_C'  ,'Inertia_Continuous'                              ,{'pointMass','plot'};...
                'sM_C'  ,'MassOffset_Continuous'                           ,{'pointMass','plot'};...
                'sMPS_C','MassPerSpan_Continuous'                          ,{'all','plot'};...
                'sMPS_C','InertiaPerSpan_Continuous'                       ,{'all','plot'};...
                'sMPS_C','MassPerSpanOffset_Continuous'                    ,{'all','plot'};...
                'sM_D'  ,'Mass_Discrete'                                   ,{'pointMass','discrete','plot'};...
                'sM_D'  ,'Inertia_Discrete'                                ,{'pointMass','discrete','plot'};...
                'sM_D'  ,'MassOffset_Discrete'                             ,{'pointMass','discrete','plot'}};
        end
        
        function NBS_Part_obj = NBS_export(obj,NBS_Master,partName,varargin)
            %assume for now that this Beam Properties Object describes a flexible nonlinear part
            NBS_Part_obj = NBS_flexPart_nonlinear(NBS_Master,partName,varargin{:});
            
            obj.is_equal(obj.sBT,obj.sK,obj.sC);
            %obj.is_equal(obj.s_alpha,obj.s_chord,obj.s_bcntr);
            obj.is_equal(obj.s_chord,obj.s_bcntr);
            assert(isempty(obj.sM_C) && isempty(obj.sMPS_C),'cannot export NBS part if sM_C and sMPS_C entries are populated. Try obj.aggregate_mass_sources to remove these.')
            
            NBS_Part_obj.s = obj.sK;
            %NBS_Part_obj.s_aero = obj.s_chord;
            NBS_Part_obj.s_ms = obj.sM_D;
            NBS_Part_obj.w = (obj.s_max-obj.s_min)*0.04; %default values for w and h
            NBS_Part_obj.h = (obj.s_max-obj.s_min)*0.005; %w and h used only for plotting
            
            NBS_Part_obj.StiffnessMatrix = obj.StiffnessMatrix;
            NBS_Part_obj.DampingMatrix = obj.DampingMatrix;
            
            NBS_Part_obj.ms = obj.Mass_Discrete;
            NBS_Part_obj.massOffset_I = obj.MassOffset_Discrete;
            NBS_Part_obj.I_varTheta = obj.Inertia_Discrete;
            
            permuteVec_ = [1 2 3 1 2 3 1 2 3];
            permuteVec = permuteVec_([2,3,4]+obj.tangent_direction);
            NBS_Part_obj.E0_W = obj.BeamTriad_reference(:,permuteVec,:);
            
            NBS_Part_obj.c = obj.chord;
            NBS_Part_obj.beam_cntr = obj.beam_cntr;
            
        end
        
    end
    
    
    
    
    methods %data manipulation methods
        function change_samplePoints(obj,sampleName,resampleVec,varargin)
            %argument sampleVecName can be:
            %    the name of an existing sample vector e.g. 'sK' (all properties based off this sample vector will be resampled)
            %    the name of a property to resample e.g. 'chord'
            %    or a group keyword e.g. 'structural'. See obj.sampleVec_property_tags to details of groups
            
            
            scalingVector = get_option(varargin,'scalingVector',1); scalingVector = reshape(scalingVector,1,1,[]);
            sampleOption = get_option(varargin,'sampleOption',''); %{'split' or ''}
            
            sVec_prop_tags = obj.sampleVec_property_tags;
            
            if nargin == 1
                for i_ = 1:size(sVec_prop_tags,1)
                    %TODO replace with sprintf statement
                    disp([sVec_prop_tags(i_,1:2) sVec_prop_tags{i_,3}])
                end
                return;
            end
            
            obj_pre_sample = copy(obj);
            
            for i_ = 1:size(sVec_prop_tags,1)
                sVecName_i = sVec_prop_tags{i_,1};
                propertyName_i = sVec_prop_tags{i_,2};
                if ~isempty(obj.(propertyName_i)) && ismember(sampleName , [{sVec_prop_tags{i_,1}} sVec_prop_tags{i_,3}])%    isequal(sampleName,sVecName_i) || isequal(sampleName,'all') )
                    if ismember('discrete',sVec_prop_tags{i_,3}) || strcmpi(sampleOption,'split')
                        obj.(propertyName_i) = sampleMassSplit(obj_pre_sample.(propertyName_i),obj_pre_sample.(sVecName_i),resampleVec,'dim',3).*scalingVector;
                    else
                        obj.(propertyName_i) = interp1dim(obj_pre_sample.(propertyName_i),obj_pre_sample.(sVecName_i),resampleVec,3,'method','linear','extrapolation','extrap').*scalingVector;
                    end
                    obj.(sVecName_i) = resampleVec;
                end
            end
        end
        
        function unify_mass_distributions(obj)
            %redefine the MPS_C, M_C and M_D distributions over the same
            %common set of sample points
            
            s_base = obj.sMPS_C; %sample points of the mass per span distribution
            s_point = [obj.sM_C;obj.sM_D]; %sample points of the point mass disctributions
            [~,index] = closestValue(s_base,s_point); %find the indices of the sMPS_C sample points closest to the sM_C points
            s_base(index) = s_point; %move these closest sMPS_C points to coincide with the sM_C/sM_D points
            s_base = unique([s_base(:) ; s_point]); %<- sample vector for the aggregated mass distribution
            
            %sample all three mass sources over the same s_base distribution
            obj.change_samplePoints('sMPS_C',s_base);
            obj.change_samplePoints('sM_C',s_base,'sampleOption','split');
            obj.change_samplePoints('sM_D',s_base,'sampleOption','split');
        end
        
        function aggregate_mass_sources(obj)
            %Combine the MPS_C, M_C and M_D sampled over the same common
            %set into a single discrete mass distribution
            
            assert(isequal(obj.sMPS_C,obj.sM_C) && isequal(obj.sM_C,obj.sM_D),'Mass Per Span / Mass Continuous / Mass Discrete disctributions must all share the same sample points before invoking this method.')
            
            diff_s_base = diff(obj.sM_C);
            s_ranges = reshape(...
               ([0 ; diff_s_base(:)]+[diff_s_base(:) ; 0])/2 , 1,1,[]);
            
            Mass_sets = {...
                obj.Mass_Continuous;...
                obj.MassPerSpan_Continuous.*s_ranges;...
                obj.Mass_Discrete};
            Offset_sets = {...
                obj.MassOffset_Continuous;...
                obj.MassPerSpanOffset_Continuous;...
                obj.MassOffset_Discrete};
            Inertia_sets = {...
                obj.Inertia_Continuous;...
                obj.InertiaPerSpan_Continuous.*s_ranges;...
                obj.Inertia_Discrete};
            
            Mass_aggregated = Mass_sets{1} + Mass_sets{2} + Mass_sets{3};
            Offsets_aggregated = (...
                Mass_sets{1}.*Offset_sets{1} + ...
                Mass_sets{2}.*Offset_sets{2} + ...
                Mass_sets{3}.*Offset_sets{3} ) ./ Mass_aggregated;
            Inertia_aggregated = Inertia_sets{1} + Inertia_sets{2} + Inertia_sets{3};
            
            %add the aggregated mass properties to the discrete property group
            obj.Mass_Discrete = Mass_aggregated;
            obj.MassOffset_Discrete = Offsets_aggregated;
            obj.Inertia_Discrete = Inertia_aggregated;
        
            %set Continuous Mass properties to zero
            obj.sM_C = [];
            obj.Mass_Continuous = 0;
            obj.MassOffset_Continuous = [0;0;0];
            obj.Inertia_Continuous = zeros(3);
            obj.sMPS_C = [];
            obj.MassPerSpan_Continuous = 0;
            obj.MassPerSpanOffset_Continuous = [0;0;0];
            obj.InertiaPerSpan_Continuous = zeros(3);
            
            %Create Array of Mass Values + Offsets
            
            
            % %             strrep(string,'oldstring','newstring')
% 
%             %add all mass sources together and write to discrete fields
%             sVec_prop_tags = obj.sampleVec_property_tags;
%             for i_ = 1:size(sVec_prop_tags,1)
%                 sVecName_i = sVec_prop_tags{i_,1};
%                 propertyName_i = sVec_prop_tags{i_,2};
%                 if ismember(sVecName_i,{'sMPS_C','sM_C'})
%                     propertyName_append = strrep(propertyName_i,'PerSpan','');
%                     propertyName_append = strrep(propertyName_i,'_Continuous','_Discrete');
%                     if ~isempty(obj.(propertyName_i))
%                         obj.(propertyName_append) = bsxfun(@plus, obj.(propertyName_append) , obj.(propertyName_i));
%                         obj.(propertyName_i) = [];
%                     end
%                 end
%             end
        end
        
        function truncate_carrythrough(obj)
            nK = numel(obj.sK);
            stiffness_norms = zeros(1,1,nK);
            StiffnessMatrices = obj.StiffnessMatrix;
            for i_ = 1:nK
                stiffness_norms(i_) = norm(StiffnessMatrices(:,:,i_));
            end
            
            %plot the stiffness matrix norms and highlight the outliers
            windowLength = 0.2;
            outlierBoolArray = isoutlier(stiffness_norms,'movmedian',floor(nK*windowLength));
            f = figure; ax = gca;
            plot(ax,obj.sK(:),stiffness_norms(:));
            plot(ax,obj.sK(outlierBoolArray),stiffness_norms(outlierBoolArray),'o');
        end
    end
    

    
    methods
        function val = get.s_min(obj)
            sVec_prop_tags = obj.sampleVec_property_tags;
            sVals = [];
            for i_ = 1:size(sVec_prop_tags,1)
                sVals = [sVals ; obj.(sVec_prop_tags{i_,1})(:)];
            end
            val = min(sVals);
        end
        function val = get.s_max(obj)
            sVec_prop_tags = obj.sampleVec_property_tags;
            sVals = [];
            for i_ = 1:size(sVec_prop_tags,1)
                sVals = [sVals ; obj.(sVec_prop_tags{i_,1})(:)];
            end
            val = max(sVals);
        end
    end
    
    methods %set methods for object properties
        %------------------------------------------------------------------
        function set.BeamTriad_reference(obj,val)
            obj.dim_match(val,[3 3 Inf]);
            obj.is_orthNorm(val);
            obj.BeamTriad_reference = val;
        end
        function set.tangent_direction(obj,val)
            [a,~] = ismember(val,[1 2 3]);
            assert(numel(a)==1 & all(a) , 'value must be a single member of the set {1,2,3}');
            obj.tangent_direction = val;
        end
        function set.sBT(obj,val)
            obj.is_vector(val);
            obj.sBT = reshape(val,1,1,[]);
        end
        %------------------------------------------------------------------
        function set.StiffnessMatrix(obj,val)
            obj.dim_match(val,[6 6 Inf]);
            obj.is_symmetric(val);
            obj.StiffnessMatrix = val;
        end
        function set.sK(obj,val)
            obj.is_vector(val);
            obj.sK = reshape(val,1,1,[]);
        end
        function set.StiffnessFnc(obj,val)
            assert(isa(@f,'function_handle') , 'value must be a function handle');
            obj.StiffnessFnc = val;
        end
        %------------------------------------------------------------------
        function set.DampingMatrix(obj,val)
            obj.dim_match(val,[6 6 Inf]);
            obj.is_symmetric(val);
            obj.DampingMatrix = val;
        end
        function set.sC(obj,val)
            obj.is_vector(val);
            obj.sC = reshape(val,1,1,[]);
        end
        function set.DampingFnc(obj,val)
            assert(isa(@f,'function_handle') , 'value must be a function handle');
            obj.DampingFnc = val;
        end
        %------------------------------------------------------------------
        function set.ReferenceOffset(obj,val)
            if val~=0
                obj.dim_match(val,[3 1 Inf]);
            end
            obj.ReferenceOffset = val;
        end
        function set.sRO(obj,val)
            obj.is_vector(val);
            obj.sRO = reshape(val,1,1,[]);
        end
        %------------------------------------------------------------------
        function set.Mass_Continuous(obj,val)
            obj.is_vector(val);
            obj.is_positive(val);
            obj.Mass_Continuous = reshape(val,1,1,[]);
        end
        function set.Inertia_Continuous(obj,val)
            obj.dim_match(val,[3 3 Inf]);
            obj.is_symmetric(val);
            obj.Inertia_Continuous = val;
        end
        function set.MassOffset_Continuous(obj,val)
            if val~=0
                obj.dim_match(val,[3 1 Inf]);
            end
            obj.MassOffset_Continuous = val;
        end
        function set.sM_C(obj,val)
            obj.is_vector(val);
            obj.sM_C = reshape(val,1,1,[]);
        end
        %------------------------------------------------------------------
        function set.MassPerSpan_Continuous(obj,val)
            obj.is_vector(val);
            obj.is_positive(val);
            obj.MassPerSpan_Continuous = reshape(val,1,1,[]);
        end
        function set.InertiaPerSpan_Continuous(obj,val)
            obj.dim_match(val,[3 3 Inf]);
            obj.is_symmetric(val);
            obj.InertiaPerSpan_Continuous = val;
        end
        function set.MassPerSpanOffset_Continuous(obj,val)
            if val~=0
                obj.dim_match(val,[3 1 Inf]);
            end
            obj.MassPerSpanOffset_Continuous = val;
        end
        function set.sMPS_C(obj,val)
            obj.is_vector(val);
            obj.sMPS_C = reshape(val,1,1,[]);
        end
        %------------------------------------------------------------------
        function set.Mass_Discrete(obj,val)
            obj.is_vector(val);
            obj.is_positive(val);
            obj.Mass_Discrete = reshape(val,1,1,[]);
        end
        function set.Inertia_Discrete(obj,val)
            obj.dim_match(val,[3 3 Inf]);
            obj.is_symmetric(val);
            obj.Inertia_Discrete = val;
        end
        function set.MassOffset_Discrete(obj,val)
            if val~=0
                obj.dim_match(val,[3 1 Inf]);
            end
            obj.MassOffset_Discrete = val;
        end
        function set.sM_D(obj,val)
            obj.is_vector(val);
            obj.sM_D = reshape(val,1,1,[]);
        end
        %------------------------------------------------------------------
        
        
        
        %------------------------------------------------------------------
%         function set.alpha_A(obj,val)
%             obj.is_vector(val);
%             obj.alpha_A = reshape(val,1,1,[]);
%         end
%         function set.s_alpha(obj,val)
%             obj.is_vector(val);
%             obj.s_alpha = val;
%         end
        %------------------------------------------------------------------
        function set.chord(obj,val)
            obj.is_vector(val);
            obj.chord = reshape(val,1,1,[]);
        end
        function set.s_chord(obj,val)
            obj.is_vector(val);
            obj.s_chord = val;
        end
        %------------------------------------------------------------------
        function set.beam_cntr(obj,val)
            obj.is_vector(val);
            obj.beam_cntr = reshape(val,1,1,[]);
        end
        function set.s_bcntr(obj,val)
            obj.is_vector(val);
            obj.s_bcntr = val;
        end
        %------------------------------------------------------------------
    end
    
    methods (Static) %assert methods performing checks for specific data properties
        function vecDim = is_vector(val)
            %true if no more than 1 dimension has length greater than 1
            bool = sum(size(val)>1)<2;
            assert(bool,'value must be a vector');
            [~,vecDim] = max(size(val));
        end
        function is_symmetric(val)
            bool = isequal( val-permute(val,[2 1 3]) , val*0);
            assert(bool,'value must be symmetric matrix in the first two dimensions');
        end
        function is_positive(val)
            bool = all(val(:)>=0);
            assert(bool,'all components of value must be positive');
        end
        function dim_match(val,sizeVec)
            %compares dimensions of val to sizeVec
            %will ignore sizeVec components set to Inf
            sz = [size(val,1),size(val,2),size(val,3)];
            temp1 = (sz == sizeVec);
            temp2 = temp1(sizeVec~=Inf);
            bool = all(temp2);
            assert(bool,'invalid dimension');
        end
        function is_orthNorm(val)
            vz = size(val,3);
            boolArray = zeros(vz,1);
            tolerance = 1e-5;
            eye3 = eye(3);
            for i_ = 1:size(val,3)
                slice = val(:,:,i_);
                temp = slice.'*slice;
                boolArray(i_) = norm(temp-eye3) < tolerance;
            end
            assert( all(boolArray) ,'matrix must represent an orthonormal triad');
        end
        function is_equal(varargin)
            varargin_1 = varargin{1};
            boolSum = 0;
            for i_ = 2:numel(varargin)
                boolSum = boolSum + isequal(varargin_1,varargin{i_});
            end
            assert(boolSum == numel(varargin)-1 , 'input quantities must be equal in size and value')
        end
    end
    
    methods %plotting functions
        
        function obj = plotBeamProperties(obj)
            figure; hold on;
            sVec_prop_tags = obj.sampleVec_property_tags;
            plotBool_cell = cellfun(@(x) ismember('plot',x),sVec_prop_tags(:,3),'UniformOutput',false);
            sVec_prop_tags_plot = sVec_prop_tags(cell2mat(plotBool_cell),:);
            nplots = size(sVec_prop_tags_plot,1);
            sp_y = 3; sp_x = ceil(nplots/sp_y);
            for sp_i = 1:nplots
                xy = sVec_prop_tags_plot(sp_i,:);
                subplot(sp_x,sp_y,sp_i); ax(sp_i) = gca; hold on;
                
                ax_i = ax(sp_i);
                x = obj.(xy{1});
                Y = obj.(xy{2});
                if isempty(Y), x = []; Y = []; end
                legendCell = {};
                for i_ = 1:size(Y,1)
                    if size(Y,2)>=i_
                        j_ = i_;
                    else
                        j_ = 1;
                    end
                    y = Y(i_,j_,:);
                    if ~isempty(y) && max(abs(y))/min(abs(y))>1e2 && min(y)>=0
                        set(ax_i,'YScale','log');
                    end
                    plot(ax_i,x(:),y(:),'.-');
                    legendCell{i_} = ['(' num2str(i_) ',' num2str(j_) ')'];
                end
                xlabel(ax_i,'s(m)','Interpreter','none'); title(ax_i,xy{2},'Interpreter','none');
                legend(ax_i,legendCell,'Interpreter','none');
            end
%             %..........................................................................
%             ax_i = 1;
%             x = obj.sK;
%             y11 = obj.StiffnessMatrix(1,1,:); plot(ax(ax_i),x(:),y11(:),'.-');
%             y22 = obj.StiffnessMatrix(2,2,:); plot(ax(ax_i),x(:),y22(:),'.-');
%             y33 = obj.StiffnessMatrix(3,3,:); plot(ax(ax_i),x(:),y33(:),'.-');
%             y44 = obj.StiffnessMatrix(4,4,:); plot(ax(ax_i),x(:),y44(:),'.-');
%             y55 = obj.StiffnessMatrix(5,5,:); plot(ax(ax_i),x(:),y55(:),'.-');
%             y66 = obj.StiffnessMatrix(6,6,:); plot(ax(ax_i),x(:),y66(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Kii');
%             legend(ax(ax_i),'K11','K22','K33','K44','K55','K66');
%             %..........................................................................
%             ax_i = 3;
%             x = obj.sM_C;
%             y = obj.Mass_Continuous; plot(ax(ax_i),x(:),y(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Mass Per Span');
%             %..........................................................................
%             ax_i = 4;
%             x = obj.sMPS;
%             y11 = obj.InertiaPerSpan(1,1,:); plot(ax(ax_i),x(:),y11(:),'.-');
%             y22 = obj.InertiaPerSpan(2,2,:); plot(ax(ax_i),x(:),y22(:),'.-');
%             y33 = obj.InertiaPerSpan(3,3,:); plot(ax(ax_i),x(:),y33(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Inertia Per Span');
%             legend(ax(ax_i),'I11','I22','I33');
%             %..........................................................................
%             ax_i = 3;
%             x = obj.sMPS;
%             y = obj.MassPerSpan; plot(ax(ax_i),x(:),y(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Mass Per Span');
%             %..........................................................................
%             ax_i = 4;
%             x = obj.sMPS;
%             y11 = obj.InertiaPerSpan(1,1,:); plot(ax(ax_i),x(:),y11(:),'.-');
%             y22 = obj.InertiaPerSpan(2,2,:); plot(ax(ax_i),x(:),y22(:),'.-');
%             y33 = obj.InertiaPerSpan(3,3,:); plot(ax(ax_i),x(:),y33(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Inertia Per Span');
%             legend(ax(ax_i),'I11','I22','I33');
%             %..........................................................................
%             ax_i = 5;
%             x = obj.sMD;
%             y = obj.MassDiscrete; plot(ax(ax_i),x(:),y(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Mass Discrete');
%             %..........................................................................
%             ax_i = 6;
%             x = obj.sMD;
%             y11 = obj.InertiaDiscrete(1,1,:); plot(ax(ax_i),x(:),y11(:),'.-');
%             y22 = obj.InertiaDiscrete(2,2,:); plot(ax(ax_i),x(:),y22(:),'.-');
%             y33 = obj.InertiaDiscrete(3,3,:); plot(ax(ax_i),x(:),y33(:),'.-');
%             set(ax(ax_i),'YScale','log');
%             xlabel(ax(ax_i),'s(m)'); title(ax(ax_i),'Inertia Discrete');
%             legend(ax(ax_i),'I11','I22','I33');
            
        end
        
    end
    
end








%     methods (Static) %data manipulation methods, static operations
%         function data_new = resample(data,s_old,s_new,dim,varargin)
%             %resample data along the indicated dimension 'dim'
%             
%             method = get_option(varargin,'method','linear'); %choose from 'nearest', 'next', 'previous', 'linear','spline','pchip', 'makima', or 'cubic'
%             extrapolation = get_option(varargin,'extrapolation','extrap'); %choose 'extrap' or provide a scalar value
%             
%             permuteVec_temp = [1 2 3];
%             permuteVec_temp(dim) = [];
%             permuteVec = [dim permuteVec_temp];
%             
%             data_permute = permute(data,permuteVec);
%             
%             data_new_permute = interp1(squeeze(s_old),data_permute,squeeze(s_new),'method',method,'extrapolation',extrapolation);
%             
%             data_new = permute(data_new_permute,permuteVec);
%         end
%     end
