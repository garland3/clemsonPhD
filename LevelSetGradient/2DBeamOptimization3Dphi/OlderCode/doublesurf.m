function h = doublesurf(varargin)
%DOUBLESURF: Produce 2 surfaces with different colormaps
%
%doublesurf(x1,y1,z1,cm1, x2, y2, z2, cm2)
% Produces 2 surfaces with different colormaps.
% x1, y1, z1 define surface 1, and is coloured with map cm1 (Nx3)
% x2, y2, z2 define surface 2, which is coloured with map cm2 (Mx3)
%
%doublesurf(z1,cm1, z2, cm2)
% Behaves as expected (z1, z2, mapped to indices).

% Copyright (c) 2001 by OPTI-NUM solutions

switch nargin
case 4
    [z1, cm1, z2, cm2] = deal(varargin{:});
    x1=1:size(z1,1);
    y1=1:size(z1,2);
    x2=1:size(z2,1);
    y2=1:size(z2,2);
case 8
    [x1, y1, z1, cm1, x2, y2, z2, cm2] = deal(varargin{:});
otherwise
    error('Incorrect calling syntax. I''m pretty sticky on this!');
end

if size(cm1,2)~=3 | size(cm2,2)~=3
    error('Colormaps must have 3 columns.');
end
% Now we make up the color indices.
zfloor1=linspace(min(z1(:)),max(z1(:)), size(cm1,1));
cind1=zeros(size(z1));
for k=1:prod(size(z1)), cind1(k)=sum(zfloor1<=z1(k));end

% The second indices must be offset by the size of the first colormap
zfloor2=linspace(min(z2(:)),max(z2(:)), size(cm2,1));
cind2=zeros(size(z2));
for k=1:prod(size(z2)), cind2(k)=sum(zfloor2<=z2(k))+size(cm1,1);end

% And the actual colormap
cmap = [cm1;cm2];

% And now we make the surfaces:
holdState = ishold;
h1 = surf(x1, y1, z1, cind1, 'CDataMapping','direct');
hold on;
h2 = surf(x2, y2, z2, cind2, 'CDataMapping','direct','edgecolor','none');
if ~holdState,
    hold off;
end
colormap(cmap);

if nargout,
    h = [h1 h2];
end