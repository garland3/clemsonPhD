
function [gauss_points, gauss_weights]=GetQuadGauss3x3()
  
  gauss_points=zeros(8,3);
  gauss_weights=zeros(8,1);
  %---first point
  gauss_points(1,1)=1.0/sqrt(3);
  gauss_points(1,2)=1.0/sqrt(3);
  gauss_points(1,3)=1.0/sqrt(3);
  gauss_weights(1)=1.0;
  
  %---second point  
  gauss_points(2,1)=-1.0/sqrt(3);
  gauss_points(2,2)=1.0/sqrt(3);
  gauss_points(2,3)=1.0/sqrt(3);
  gauss_weights(2)=1.0;
  
  %---third point  
  gauss_points(3,1)=-1.0/sqrt(3);
  gauss_points(3,2)=-1.0/sqrt(3);
  gauss_points(3,3)=1.0/sqrt(3);
  gauss_weights(3)=1.0;
 
  %---fourth point  
  gauss_points(4,1)=1.0/sqrt(3);
  gauss_points(4,2)=-1.0/sqrt(3);
  gauss_points(4,3)=1.0/sqrt(3);
  gauss_weights(4)=1.0;
  
  %---fifth point
  gauss_points(5,1)=1.0/sqrt(3);
  gauss_points(5,2)=1.0/sqrt(3);
  gauss_points(5,3)=-1.0/sqrt(3);
  gauss_weights(5)=1.0;
  
  %---sixth point  
  gauss_points(6,1)=-1.0/sqrt(3);
  gauss_points(6,2)=1.0/sqrt(3);
  gauss_points(6,3)=-1.0/sqrt(3);
  gauss_weights(6)=1.0;
  
  %---seventh point  
  gauss_points(7,1)=-1.0/sqrt(3);
  gauss_points(7,2)=-1.0/sqrt(3);
  gauss_points(7,3)=-1.0/sqrt(3);
  gauss_weights(7)=1.0;
 
  %---eighth point  
  gauss_points(8,1)=1.0/sqrt(3);
  gauss_points(8,2)=-1.0/sqrt(3);
  gauss_points(8,3)=-1.0/sqrt(3);
  gauss_weights(8)=1.0;
      