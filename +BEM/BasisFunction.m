classdef BasisFunction < BEUT.BEM.PiecewisePolynomial
    % Define the basis function as a piecewise polynomial on every segment of the geometry
    
    properties
        pol;          % polynomial basis function(s) for each edge
        idx;          % polynomial indices to define which polynomial segments to apply
        
        idx_table;    % index table to define which polynomial segments to apply to each edge 
    end
    
    methods
        % Constructor
        function obj = BasisFunction()
            
        end % constructor
        
    end % methods
    
    methods(Static)
        
        % Create hat basis, "scale" set to true scales by 1/length
        function obj = createHat(halfedges, scale)
            
            if nargin<2, scale=false; end
            degree = 1;
            
            N_E = length(halfedges);
            shape_num = vertcat(halfedges.shape);
            num_shapes = max(shape_num);
            
            % 2 polynomial contributions per segment
            N_S = 2;
            obj.pol = cell(N_S,N_E);
            
            for edge = 1 : N_E
                
                if scale==true
                    l = halfedges(edge).l;
                else l=1;
                end
                
                % upslope = s
                obj.pol{1,edge} = BEUT.BEM.PiecewisePolynomial([0 1], [1 0]/l, degree);
                % downslope = -s+1
                obj.pol{2,edge} = BEUT.BEM.PiecewisePolynomial([0 1], [-1 1]/l, degree);
                
                obj.idx{edge} = [1 2];
            end
            
            % Because there may be more than 1 shape, we don't want the basis functions connecting spatially
            % distinct shapes, i.e. they need to loop around each individual object
            shape_N_E = 0;
            for shape=1:num_shapes
                shape_N_E = [shape_N_E numel(find(shape_num==shape))];
            end
            shape_N_E(1) = [];
            
            obj.idx_table = create_index_table(shape_N_E,N_S,shape_N_E,1);
            
        end % createHat
        
        % Create square basis, "scale" set to true scales by 1/length
        function obj = createSquare(halfedges,scale)
            
            if nargin<2, scale=true; end
            degree = 0;
            
            numEdges = length(halfedges);
            
            % 1 polynomial contribution per segment
            obj.pol = cell(1,numEdges);
            
            for edge = 1 : numEdges
                
                if scale==true
                    l = halfedges(edge).l;
                else l=1;
                end
                
                % constant
                obj.pol{1,edge} = BEUT.BEM.PiecewisePolynomial([0 1], 1/l, degree);
                
                obj.idx{edge} = 1;
            end
            
            obj.idx_table = create_index_table(numEdges,1,1,1);
            
        end % createSquare
        
        % Create dual square basis, "scale" set to true scales by 1/length
        function obj = createDualSquare(dual_halfedges, scale)
            
            if nargin<2, scale=true; end
            degree = 0;
            
            N_E = length(dual_halfedges);
            assert(~logical(bitand(N_E,1)),['Dual mesh will always have an even number of edges, '...
                'the geometry input into this function has an odd number of edges.']);
            
            % 2 polynomial contribution per function
            N_S = 2;
            obj.pol = cell(N_S,N_E);
            
            % zero polynomial
            zero_pol = BEUT.BEM.PiecewisePolynomial([0,1], 0, 0);
            
            for edge = 2 : 2 : N_E
                
                % Each edge is split into 2, thus the pulse will need to be doubled
                % to get back to original length
                if scale==true
                    l = dual_halfedges(edge).l * 2;
                else l=1;
                end
                
                % constant
                obj.pol{1,edge-1} = BEUT.BEM.PiecewisePolynomial([0 1], 1/l, degree);
                obj.pol{2,edge-1} = zero_pol;
                obj.pol{1,edge} = zero_pol;
                obj.pol{2,edge} = BEUT.BEM.PiecewisePolynomial([0 1], 1/l, degree);
                
                obj.idx{edge-1} = 1;
                obj.idx{edge} = 2;
            end
            
            obj.idx_table = create_index_table(N_E,N_S,1,2);
            
        end % createDualSquare
        
        
        function obj = createDualHat(dual_halfedges, scale)
            
            if nargin<2, scale=false; end
            degree = 1;
            
            N_E = length(dual_halfedges);
            shape_num = vertcat(dual_halfedges.shape);
            num_shapes = max(shape_num);
            
            assert(~logical(bitand(N_E,1)),['Dual mesh will always have an even number of edges, '...
                'the geometry input into this function has an odd number of edges.']);
            
            % 4 polynomial contributions per function
            N_S = 4;
            obj.pol = cell(N_S,N_E);
            
            % zero polynomial
            zero_pol = BEUT.BEM.PiecewisePolynomial([0,1], [0 0], 0);
            
            for edge = 2 : 2 : N_E
                
                % Each edge is split into 2, thus the pulse will need to be doubled
                % to get back to original length
                if scale==true
                    l = dual_halfedges(edge).l * 2;
                else l=1;
                end
                
                % upslope 1 = s/2
                obj.pol{1,edge} = BEUT.BEM.PiecewisePolynomial([0 1], [0.5 0]/l, degree);
                obj.pol{1,edge-1} = zero_pol;
                % upslope 2 = s/2 + 1/2
                obj.pol{2,edge} = zero_pol;
                obj.pol{2,edge-1} = BEUT.BEM.PiecewisePolynomial([0 1], [0.5 0.5]/l, degree);
                % downslope 3 = -s/2 + 1
                obj.pol{3,edge} = BEUT.BEM.PiecewisePolynomial([0 1], [-0.5 1]/l, degree);
                obj.pol{3,edge-1} = zero_pol;
                % downslope 4 = -s/2 + 1/2
                obj.pol{4,edge} = zero_pol;
                obj.pol{4,edge-1} = BEUT.BEM.PiecewisePolynomial([0 1], [-0.5 0.5]/l, degree);        
                
                obj.idx{edge} = [1 3];
                obj.idx{edge-1} = [2 4];
                
            end
            
            % Because there may be more than 1 shape, we don't want the basis functions connecting spatially
            % distinct shapes, i.e. they need to loop around each individual object
            shape_N_E = 0;
            for shape=1:num_shapes
                shape_N_E = [shape_N_E numel(find(shape_num==shape))];
            end
            shape_N_E(1) = [];
            
            obj.idx_table = create_index_table(shape_N_E,N_S,shape_N_E,2);
            
        end % createDualHat
        
        
        % Divergence operator
        function obj = divergence(obj, halfedges)
            
            numBasisFunctions = size(obj.pol,2);
            assert(numBasisFunctions==numel(halfedges))
            
            numElements = size(obj.pol,1);
            
            for segment = 1:numBasisFunctions
                
                l = halfedges(segment).l;
                
                for e = 1:numElements
                    
                    obj.pol{e,segment} = diff(obj.pol{e,segment});

                    obj.pol{e,segment}.coeffs = obj.pol{e,segment}.coeffs / l;
                    
                end
                
                
            end
            
        end % div
        
        
        % Collapse the 2D geometry into a horizontal line and plot
        % basis functions associated with each edge on that line
        function plot_basis(halfedges,basis)
            
            basis = vertcat(basis.pol);
            
            N_V = size(basis,2);
            num_basis = size(basis,1);
            
            assert(N_V==length(halfedges),['Both inputs to plot_basis must relate to the same geometry,'...
                'i.e. have the same number of edges']);
            
            figure; hold on;
            axis equal;
            % color array for cycling through basis functions
            color = [   1 0 0;...      red
                        0 0 1;...      blue
                        0 1 0;...      green
                        1 0 1;...      magenta
                        0 1 1;...      cyan
                        0.6 0.4 0.2;...brown
                        0 0.5 0;...    dark green
                        1 0.5 0;...    orange
                        0 0 0;...      black
                        0.9 0.9 0;...  yellow
                        0.5 0 0.5;]; % purple
            
            % accuracy of each basis
            x_ax = 0:0.01:1;
            
            start_edge = 0;
            for edge_num=1:N_V
                l = halfedges(edge_num).l;
                
                for basis_num=1:num_basis
                    
                    % Geometry is parameterised as a straight horizontal line
                    plot([start_edge,start_edge+l],[0,0],'o-k', 'LineWidth',2)
                    
                    % Basis plot
                    plot(x_ax*l+start_edge, basis{basis_num,edge_num}.eval(x_ax),...
                        'Color',color(basis_num,:));
                
                end
                start_edge = start_edge + l;
            end
            
        end % plot
        
        
    end % static methods
    
end

% Create index table to decide the basis polynomial segments that apply to each edge
function table = create_index_table(N_E, N_S, start, iterate)
% N_E = number of edges (if more than 1 shape, this can be an array, with the sum being total number of edges)
% N_S = number of function segments
% start = edge number of which the first function starts (can be an array if more than one shape)
% iterate = how many edges to skip to iterate to the next function, i.e. if iterate=1 then a function will
%           begin at each edge, iterate=2 and a function will begin every 2nd edge etc.

assert(numel(N_E)==numel(start), 'Array for N_E and start must be consistent');
num_shapes = numel(N_E);

% find nunmber of shapes
if num_shapes==1
    
    table = create_index_table_for_each_shape(N_E,N_S,start,iterate);
    
else
    
    % first shape
    table = create_index_table_for_each_shape(N_E(1),N_S,start(1),iterate);
    
    % all other shapes
    for j=2:num_shapes
        
        tbl_temp = create_index_table_for_each_shape(N_E(j),N_S,start(j),iterate);
        
        table = [table; tbl_temp+N_E(j-1)];
        
    end
    
end

    function table = create_index_table_for_each_shape(N_E,N_S,start,iterate)
        a = (1:N_E)'-2 + start;
        
        mod_=@(array,N) mod(array-1,N)+1;
        
        b = zeros(N_E,N_S);
        for i=1:N_S
            b(:,i)=mod_(a+i,N_E);
        end
        
        filter = 1:iterate:N_E;
        table = b(filter,:);
    end
end

