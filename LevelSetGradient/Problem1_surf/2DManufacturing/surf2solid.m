function varargout = surf2solid(varargin)
%SURF2SOLID  Convert a thin surface to a closed triangulated solid volume.
%
%   SOLID_FV = SURF2SOLID(FV,...) takes in a triangulated patch defined by
%   FV (a structure with fields 'vertices' and 'faces'), and returns a
%   solid patch SOLID_FV closed by options (described below).
%
%   SOLID_FV = SURF2SOLID(F, V,...) takes faces and vertices separately.
%
%   [F,V] = SURF2SOLID(...) returns solid faces and vertices separately.
%
%   SURF2SOLID(...) with no output argument plots the 3 components
%   (orig-surface, side-walls, under-surface) to a new figure.
%
%   SURF2SOLID(X, Y, Z, ...) reads in surface data in X, Y, and Z matrices,
%   and triangulates this gridded data into a surface using triangulation
%   options specified below. Z must be a 2D matrix. X and Y can be 2D
%   matrices of the same size as Z, or vectors of length equal to SIZE(Z,2)
%   and SIZE(Z,1), respectively. If X or Y are scalar values, they are used
%   to specify the X and Y spacing between grid points.
%
%   SURF2SOLID(...,'PropertyName',VALUE,...) makes a solid volume from thin
%   surface using any of the following property/value options:
%
%   ELEVATION     - Extends the surface down to a flat base at the given
%                   (Z) elevation value. Useful for turning a thin
%                   elevation map into a solid block with a flat base. The
%                   ELEVATION value should be below the lowest (or above
%                   the highest) data point. If no other options are given,
%                   ELEVATION defaults to MIN(Z)-0.1*(MAX(Z)-MIN(Z)).
%                   Variable ELEVATION may also be given per-point, via a
%                   2D matrix (the same size as Z for X,Y,Z style input) or
%                   a 1D array (with length equal to the number of vertices
%                   given in face/vertex input).
%
%   THICKNESS     - Value to offset the given thin surface to make a
%                   thickened solid slab. Each node on the surface will be
%                   projected along its normal direction by thickness. When
%                   negative thickness is given, offset will be away from
%                   face normal direction. Variable thickness can also be
%                   specified via a 2D matrix (of same size as Z, for X,Y,Z
%                   input) or an N-by-1 array of thicknesses (where N is
%                   the number of vertices in the thin surface)
%
%   TRIANGULATION - When used with gridded data, TRIANGULATION is either:
%                    'delaunay'  - (default) Delaunay triangulation of X, Y
%                    'f'         - Forward slash division of grid quads
%                    'b'         - Back slash division of quadrilaterals
%                    'x'         - Cross division of quadrilaterals
%                   Note that 'f', 'b', or 'x' triangulations use an
%                   inbuilt version of FEX entry 28327, "mesh2tri". 'x'
%                   style triangulation cannot be used with variable
%                   ELEVATION or THICKNESS parameters.
%
%   NORMALS       - When THICKNESS options is used, the direction to
%                   thicken the surface is (by default) determined by the
%                   surface (unit vector) normal directions at each vertex.
%                   To override these default directions, you may specify
%                   NORMALS as an N-by-3 array of normal directions (where
%                   N is the number of vertices in the thin surface). This
%                   is useful when underlying data gives more precise
%                   normal directions than face orienatations (for an
%                   example, see the isonormals function).
%
% Note 1: Currently surf2solid will return a closed surface with face
% normals pointing "out". With user feedback, I'd be happy to change this
% behaviour to either "in" or "unchanged from input direction".
% Note 2: If a single ELEVATION value is specified (i.e., flat base), the
% resulting patch will have minimal triangles on the flat base to reduce
% patch/file size.
%
%   Example (shows both THICKNESS and ELEVATION forms):
%     n = 30;
%     [X,Y] = meshgrid(linspace(0,1,2*n+1));
%     L = (40/51/0.9)*membrane(1,n);
%     figure, subplot(2,2,[1 3]), title 'Thin surface'
%     surf(X,Y,L,'EdgeColor','none'); colormap pink; axis image; camlight
%     subplot(2,2,2), title 'Block elevation'
%     surf2solid(X,Y,L,'elevation',min(L(:))-0.05); axis image; camlight; camlight
%     subplot(2,2,4), title 'Thickness'
%     surf2solid(X,Y,L,'thickness',-0.1); axis image; camlight;
%
%   Original idea adapted from Paul Kassebaum's blog post
%   http://blogs.mathworks.com/community/2013/06/20/paul-prints-the-l-shaped-membrane/
%   Many thanks to Paul for his further input and improvements.
%
%   Author: Sven Holcombe, 07-20-2013

% 1.0 (2013-07) Original
% 1.01(2013-07) Set other-faces as opposite norm direction to input faces
%               Ensured optional per-node thickness input is allowed
% 1.1 (2013-07) Added per-node elevation and thickness, minimal-flat-base
%               file size, and default face orientation.
% 1.2 (2013-08) Reduced-mode limited to gridded input.
% 1.3 (2014-02) Added optional face normals input.

% Get faces, vertices, and user-defined options for writing
[F, V, options] = parseInputs(varargin{:});

% Get the latest triangulation class. 2013a brought in "triangulation" with
% some nice options. Earlier versions must use "TriRep".
if exist('triangulation','class')
    T = triangulation(F,V);
    options.oldVersion = false;
else
    T = TriRep(F,V); %#ok<DTRIREP>
    options.oldVersion = true;
end
% Extract boundary edges from input surface. These will connect to walls.
boundEdges = T.freeBoundary;
if ~isempty(boundEdges)
    boundVerts = boundEdges([1:end 1],1);
else
    boundVerts = zeros(0,1);
end

% Define "other" faces opposite to input faces. These will either sit at
% given Z-elevations, or they will be offset from input faces by thickness.
if ~isempty(options.elevation)  % ELEVATION was specified
    % If scalar elevation was given, it will be assigned here. If variable
    % elevation was given, it will also be assigned.
    V_extrude = V;
    V_extrude(:,3) = options.elevation(:);
else                            % THICKNESS was specified
    % Vertex normals may be provided as input ...
    if ~isempty(options.normals)
        Vnormals = options.normals;
        % Or we calculate vertex normals the hard way ...
    elseif options.oldVersion
        facets = V';
        facets = permute(reshape(facets(:,F'), 3, 3, []),[2 1 3]);
        allEdgeVecs = facets([2 3 1],:,:) - facets(:,:,:);
        allFacetNormals =  bsxfun(@times, allEdgeVecs(1,[2 3 1],:), allEdgeVecs(2,[3 1 2],:)) - ...
            bsxfun(@times, allEdgeVecs(2,[2 3 1],:), allEdgeVecs(1,[3 1 2],:));
        allFacetNormals = bsxfun(@rdivide, allFacetNormals, sqrt(sum(allFacetNormals.^2,2)));
        facesByVertex = T.vertexAttachments;
        Vnormals = zeros(size(V));
        for i = 1:length(facesByVertex)
            Vnormals(i,:) =  mean(allFacetNormals(:,:,facesByVertex{i}),3);
        end
        Vnormals = bsxfun(@rdivide, Vnormals, sqrt(sum(Vnormals.^2,2)));
    else
        % Or we calculate vertex normals the easy way
        Vnormals = T.vertexNormal;
    end
    % Extrudes by thickness in each normal direction.
    % bsxfun is used in case the user wants to supply variables offsets by
    % each vertex.
    V_extrude = V + bsxfun(@times, Vnormals, options.thickness(:));
end

% If a scalar elevation was supplied, we can save file size (or triangle
% count) and define the base by only its boundary vertices. Note that if
% faces/vertices input was supplied, it's possible a complex geometry was
% given that won't suit our reduction method below. Let's be conservative.
if isscalar(options.elevation) && options.griddedInput
    options.baseReduced = true;
    V_wall = [V(boundVerts,:); V_extrude(boundVerts,:)];
else % thickness was specified or complex geometry was possible
    options.baseReduced = false;
    V_wall = [V(boundVerts,:); V_extrude(boundVerts,:)];
end

% Number of wall vertices on each surface (nwv).
nwv = length(V_wall)/2;
% Allocate memory for wallFaces.
F_wall = zeros(2*(nwv-1),3);
% Define the faces.
for k = 1:nwv-1
    F_wall(k      ,:) = [k+1  ,k      ,k+nwv];
    F_wall(k+nwv-1,:) = [k+nwv,k+1+nwv,k+1  ];
end

% Let's use the first vertex to test if faces are pointed in/out
testNormal = cross(...
    V(F(1,2),:) - V(F(1,1),:),...
    V(F(1,3),:) - V(F(1,1),:));
if ~isempty(options.elevation)
    firstVertZoffset = options.elevation(1) - V(1,3);
else
    firstVertZoffset = options.thickness(1);
end
% If the first face is pointing in the same direction as the direction of
% offset, we will have all faces pointing "in". We want them pointing "out"
if sign(testNormal(3)) == sign(firstVertZoffset)
    F = fliplr(F);
end

if options.baseReduced       % SCALAR ELEVATION was specified
    % Each boundary vertex forms a triangle with its right-hand neighbor
    % and a newly formed centroid of the extruded face.
    n = length(boundVerts);
    V_extrude = [mean(V_extrude,1); V_extrude(boundVerts,:)];% prepend centroid.
    % Ensure extruded faces are properly oriented.
    testNormal = cross(...
        V_extrude(2,:)-V_extrude(1,:),...
        V_extrude(3,:)-V_extrude(1,:));
    if sign(testNormal(3)) == sign(options.elevation)
        F_extrude = [ones(n-1,1),(2:n)',[(3:n)';2]];
    else
        F_extrude = [[(3:n)';2],(2:n)',ones(n-1,1)];
    end
else % Thickness or variable elevation was specified
    % Simply ensure the extruded faces are oriented opposite the originals
    F_extrude = fliplr(F);
end

% Compile 3 sets of faces together
allVertices = [V; V_wall; V_extrude];
allFaces = [F;                           % Use original faces
    F_wall+size(V,1);                    % Add wall faces
    F_extrude+size(V,1)+size(V_wall,1)]; % Add opposite faces (flipped)

% Ouput based on requested variables.
varargout = {};
if nargout == 0
    %figure;
    hold on;
    view(3);
    axis vis3d;
    patch('Faces',F        , 'Vertices',V        , 'FaceColor','r');
    patch('Faces',F_wall   , 'Vertices',V_wall   , 'FaceColor','g');
    patch('Faces',F_extrude, 'Vertices',V_extrude, 'FaceColor','b');
    hold off;
    %   set(gca,'zlim',[min(allVertices(:,3)),max(allVertices(:,3))]);
elseif nargout == 1
    varargout = {struct('faces',allFaces,'vertices',allVertices)};
elseif nargout >= 2
    varargout = {allFaces, allVertices};
end


%% Input handling subfunctions
function [faces, vertices, options] = parseInputs(varargin)
% Determine input type
if isstruct(varargin{1}) % surf2solid(FVstruct, ...)
    if ~all(isfield(varargin{1},{'vertices','faces'}))
        error( 'Variable p must be a faces/vertices structure' );
    end
    faces = varargin{1}.faces;
    vertices = varargin{1}.vertices;
    options = parseOptions(varargin{2:end});
    
elseif isnumeric(varargin{1})
    firstNumInput = cellfun(@isnumeric,varargin);
    firstNumInput(find(~firstNumInput,1):end) = 0; % Only consider numerical input PRIOR to the first non-numeric
    numericInputCnt = nnz(firstNumInput);
    
    options = parseOptions(varargin{numericInputCnt+1:end});
    switch numericInputCnt
        case 3 % surf2solid(X, Y, Z, ...)
            options.griddedInput = true;
            % Extract the matrix Z
            Z = varargin{3};
            
            % Convert scalar XY to vectors
            ZsizeXY = fliplr(size(Z));
            for i = 1:2
                if isscalar(varargin{i})
                    varargin{i} = (0:ZsizeXY(i)-1) * varargin{i};
                end                    
            end
            
            % Extract X and Y
            if isequal(size(Z), size(varargin{1}), size(varargin{2}))
                % X,Y,Z were all provided as matrices
                [X,Y] = varargin{1:2};
            elseif numel(varargin{1})==ZsizeXY(1) && numel(varargin{2})==ZsizeXY(2)
                % Convert vector XY to meshgrid
                [X,Y] = meshgrid(varargin{1}, varargin{2});
            else
                error('surf2solid:badinput', 'Unable to resolve X and Y variables');
            end
            
            % Convert to faces/vertices
            if strcmp(options.triangulation,'delaunay')
                faces = delaunay(X,Y);
                vertices = [X(:) Y(:) Z(:)];
            else
                if ~exist('mesh2tri','file')
                    error('surf2solid:missing', '"mesh2tri" is required to convert X,Y,Z matrices to STL. It can be downloaded from:\n%s\n',...
                        'http://www.mathworks.com/matlabcentral/fileexchange/28327')
                end
                [faces, vertices] = mesh2tri(X, Y, Z, options.triangulation);
            end
            
        case 2 % surf2solid(FACES, VERTICES, ...)
            faces = varargin{1};
            vertices = varargin{2};
            
        otherwise
            error('surf2solid:badinput', 'Unable to resolve input types.');
    end
end
% Ensure *some* information is there to make a thickness
if isempty(options.thickness) && isempty(options.elevation)
    options.elevation = min(vertices(:,3)) - 0.1 * (max(vertices(:,3))-min(vertices(:,3)));
end

function options = parseOptions(varargin)
IP = inputParser;
IP.addParamValue('triangulation', 'delaunay', @ischar);
IP.addParamValue('elevation',[],@isnumeric)
IP.addParamValue('thickness',[],@isnumeric)
IP.addParamValue('normals',  [],@isnumeric)
IP.addParamValue('griddedInput',false)
IP.parse(varargin{:});
options = IP.Results;

%% INCLUDED FUNCTIONS %%
function [F,V]=mesh2tri(X,Y,Z,tri_type)
% function [F,V]=mesh2tri(X,Y,Z,tri_type)
% 
% Available from http://www.mathworks.com/matlabcentral/fileexchange/28327
% Included here for convenience. Many thanks to Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/07/2010
%------------------------------------------------------------------------
[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1);
switch tri_type
    case 'f'%Forward slash
        TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
        TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'b'%Back slash
        TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
        TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'x'%Cross
        TRI_I=[I(:)+1,I(:);  I(:)+1,I(:)+1;  I(:),I(:)+1;    I(:),I(:)];
        TRI_J=[J(:),J(:);    J(:)+1,J(:);    J(:)+1,J(:)+1;  J(:),J(:)+1];
        IND=((numel(X)+1):numel(X)+prod(size(X)-1))';
        F = sub2ind(size(X),TRI_I,TRI_J);
        F(:,3)=repmat(IND,[4,1]);
        Fe_I=[I(:),I(:)+1,I(:)+1,I(:)]; Fe_J=[J(:),J(:),J(:)+1,J(:)+1];
        Fe = sub2ind(size(X),Fe_I,Fe_J);
        Xe=mean(X(Fe),2); Ye=mean(Y(Fe),2);  Ze=mean(Z(Fe),2);
        X=[X(:);Xe(:)]; Y=[Y(:);Ye(:)]; Z=[Z(:);Ze(:)];
end
V=[X(:),Y(:),Z(:)];