

if(1==0)
    % Example 1: Decompose 200 random points into bins of 20 points or less,
    %             then display each bin with its points in a separate colour.
       pts = (rand(200,3)-0.5).^2;
       OT = OcTree(pts,'binCapacity',20);        
       figure
       boxH = OT.plot;
       cols = lines(OT.BinCount);
       doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
       for i = 1:OT.BinCount
           set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
           doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
       end
       axis image, view(3)
%

else
    % example 2
      pts = rand(200,3);
       OT = OcTree(pts,'binCapacity',10,'style','weighted');
       OT.shrink
       figure
       boxH = OT.plot;
       cols = lines(OT.BinCount);
       doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
       for i = 1:OT.BinCount
           set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
           doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
       end
       axis image, view(3)

end