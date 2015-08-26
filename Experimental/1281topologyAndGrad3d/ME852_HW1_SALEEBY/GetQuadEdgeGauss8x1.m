
function [gauss_points, gauss_weights]=GetQuadEdgeGauss8x1()
  
  gauss_points=zeros(24,3);
  gauss_weights=zeros(24,1);
  
  %---1st point
  gauss_points(1,1)=1.0/sqrt(3);
  gauss_points(1,2)=1.0/sqrt(3);
  gauss_points(1,3)=1;
  gauss_weights(1)=1.0;
  
  %---2nd point
  gauss_points(2,1)=-1.0/sqrt(3);
  gauss_points(2,2)=1.0/sqrt(3);
  gauss_points(2,3)=1;
  gauss_weights(2)=1.0;
  
  %---3rd point
  gauss_points(3,1)=-1.0/sqrt(3);
  gauss_points(3,2)=-1.0/sqrt(3);
  gauss_points(3,3)=1;
  gauss_weights(3)=1.0;

  %---4th point
  gauss_points(4,1)=1.0/sqrt(3);
  gauss_points(4,2)=-1.0/sqrt(3);
  gauss_points(4,3)=1;
  gauss_weights(4)=1.0;
  
  %---5th point
  gauss_points(5,1)=1.0/sqrt(3);
  gauss_points(5,2)=1.0/sqrt(3);
  gauss_points(5,3)=-1;
  gauss_weights(5)=1.0;
  
  %---6th point
  gauss_points(6,1)=-1.0/sqrt(3);
  gauss_points(6,2)=1.0/sqrt(3);
  gauss_points(6,3)=-1;
  gauss_weights(6)=1.0;
  
  %---7th point
  gauss_points(7,1)=-1.0/sqrt(3);
  gauss_points(7,2)=-1.0/sqrt(3);
  gauss_points(7,3)=-1;
  gauss_weights(7)=1.0;
  
  %---8th point
  gauss_points(8,1)=1.0/sqrt(3);
  gauss_points(8,2)=-1.0/sqrt(3);
  gauss_points(8,3)=-1;
  gauss_weights(8)=1.0;
 
  gauss_points(9,1)=1.0/sqrt(3);
  gauss_points(9,2)=1.0;
  gauss_points(9,3)=1.0/sqrt(3);
  gauss_weights(9)=1.0;
  
  gauss_points(10,1)=-1.0/sqrt(3);
  gauss_points(10,2)=1.0;
  gauss_points(10,3)=1.0/sqrt(3);
  gauss_weights(10)=1.0;
  
  gauss_points(11,1)=-1.0/sqrt(3);
  gauss_points(11,2)=1.0;
  gauss_points(11,3)=-1.0/sqrt(3);
  gauss_weights(11)=1.0;
  
  gauss_points(12,1)=1.0/sqrt(3);
  gauss_points(12,2)=1.0;
  gauss_points(12,3)=-1.0/sqrt(3);
  gauss_weights(12)=1.0;
  
  gauss_points(13,1)=1.0/sqrt(3);
  gauss_points(13,2)=-1.0;
  gauss_points(13,3)=1.0/sqrt(3);
  gauss_weights(13)=1.0;
  
  gauss_points(14,1)=-1.0/sqrt(3);
  gauss_points(14,2)=-1.0;
  gauss_points(14,3)=1.0/sqrt(3);
  gauss_weights(14)=1.0;
  
  gauss_points(15,1)=-1.0/sqrt(3);
  gauss_points(15,2)=-1.0;
  gauss_points(15,3)=-1.0/sqrt(3);
  gauss_weights(15)=1.0;
  
  gauss_points(16,1)=1.0/sqrt(3);
  gauss_points(16,2)=-1.0;
  gauss_points(16,3)=-1.0/sqrt(3);
  gauss_weights(16)=1.0;
  
  gauss_points(17,1)=1.0;
  gauss_points(17,2)=1.0/sqrt(3);
  gauss_points(17,3)=1.0/sqrt(3);
  gauss_weights(17)=1.0;
  
  gauss_points(18,1)=1.0;
  gauss_points(18,2)=-1.0/sqrt(3);
  gauss_points(18,3)=1.0/sqrt(3);
  gauss_weights(18)=1.0;
  
  gauss_points(19,1)=1.0;
  gauss_points(19,2)=-1.0/sqrt(3);
  gauss_points(19,3)=-1.0/sqrt(3);
  gauss_weights(19)=1.0;
  
  gauss_points(20,1)=1.0;
  gauss_points(20,2)=1.0/sqrt(3);
  gauss_points(20,3)=-1.0/sqrt(3);
  gauss_weights(20)=1.0;
  
  gauss_points(21,1)=-1.0;
  gauss_points(21,2)=1.0/sqrt(3);
  gauss_points(21,3)=1.0/sqrt(3);
  gauss_weights(21)=1.0;
  
  gauss_points(22,1)=-1.0;
  gauss_points(22,2)=-1.0/sqrt(3);
  gauss_points(22,3)=1.0/sqrt(3);
  gauss_weights(22)=1.0;
  
  gauss_points(23,1)=-1.0;
  gauss_points(23,2)=-1.0/sqrt(3);
  gauss_points(23,3)=-1.0/sqrt(3);
  gauss_weights(23)=1.0;
  
  gauss_points(24,1)=-1.0;
  gauss_points(24,2)=1.0/sqrt(3);
  gauss_points(24,3)=-1.0/sqrt(3);
  gauss_weights(24)=1.0;
  
  
  
  