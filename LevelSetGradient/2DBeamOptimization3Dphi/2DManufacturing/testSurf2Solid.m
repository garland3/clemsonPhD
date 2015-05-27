 n = 30;
   [X,Y] = meshgrid(linspace(0,1,2*n+1));
     L = (40/51/0.9)*membrane(1,n);
     figure, subplot(2,2,[1 3]), title 'Thin surface'
     surf(X,Y,L,'EdgeColor','none'); colormap pink; axis image; camlight
     subplot(2,2,2), title 'Block elevation'
     fv = surf2solid(X,Y,L,'elevation',min(L(:)));patch(fv); axis image; camlight; camlight
     subplot(2,2,4), title 'Thickness'
     surf2solid(X,Y,L,'thickness',-0.1); axis image; camlight;
     
     
 stlwrite('abc.stl',fv,'mode','ascii')         % Save to binary .stl