function [w] = rectangleWeight(type, para, x,y,xI,yI,dmIx,dmIy)
% EVALUATE WEIGHT FUNCTION
%
% SYNTAX:[w, dwdx, dwdy] = rectangleWeight(type, para, x,y,xI,yI,dmIx,dmIy)
%
% INPUT PARAMETERS
%    type - Type of weight function
%    para - Weight function parameter
%    x,y   - gauss point coordinates 
%    xI,yI  -  nodal point coordinate
%    dmIx - Support size with respect to x direction
%    dmIy - Support size with respect to y direction
% OUTPUT PARAMETERS
%    w    - Value of weight function at r
%    dwdx - Value of first order derivative of weight function with respect to x at r
%    dwdy - Value of first order derivative of weight function with respect to y at r
% x=repmat(x,size(xI));
% y=repmat(y,size(xI));
dmIx=reshape(dmIx,size(xI));
dmIy=reshape(dmIy,size(yI));
rx  = abs(x-xI)./dmIx;
ry  = abs(y-yI)./dmIy; %  define the support size is a rectangle

       
% EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r

if strcmp(type,'GAUSS')
   [wx] = Gauss(para,rx);
   [wy] = Gauss(para,ry);
elseif (type == 'CUBIC')
   [wx,dwdrx] = Cubic(rx);
   [wy,dwdry] = Cubic(ry);
elseif (type == 'SPLI3')
   [wx,dwdrx] = Spline3(rx);
   [wy,dwdry] = Spline3(ry);
elseif (type == 'SPLI5')
   [wx,dwdrx] = Spline5(rx);
   [wy,dwdry] = Spline5(ry);
  elseif (type == 'power')
   [wx,dwdrx] = power_function(para,rx);
   [wy,dwdry] = power_function(para,ry);
elseif (type == 'CRBF1')
   [wx,dwdrx] = CSRBF1(rx);
    [wy,dwdry] = CSRBF1(ry);
elseif (type == 'CRBF2')
   [wx,dwdrx] = CSRBF2(rx);
   [wy,dwdry] = CSRBF2(ry);
elseif (type == 'CRBF3')
   [wx,dwdrx] = CSRBF3(rx);
   [wy,dwdry] = CSRBF3(ry);
elseif (type == 'CRBF4')
   [wx,dwdrx] = CSRBF4(rx);
   [wy,dwdry] = CSRBF4(ry);
elseif (type == 'CRBF5')
   [wx,dwdrx] = CSRBF5(rx);
   [wy,dwdry] = CSRBF5(ry);
elseif (type == 'CRBF6')
   [wx,dwdrx] = CSRBF6(rx);
   [wy,dwdry] = CSRBF6(ry);
else
   error('Invalid type of weight function.');
end
   
w   =  wx.*wy;



function [w] = Gauss(beta,r)
w=r;
b2 = beta*beta;
r2 = r.*r;
eb2 = exp(-b2);
w  = (exp(-b2.*r2) - eb2) ./ (1.0 - eb2);
w(r>1.0)=0.0;



function [w,dwdr] = Cubic(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
else
   w     = 1-6*r^2+8*r^3-3*r^4;
end

function [w,dwdr] = Spline3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
elseif (r<=0.5)
   w     = 2/3 - 4*r^2 + 4*r^3;
   dwdr  = -8*r + 12*r^2;
 
else
   w     = 4/3 - 4*r + 4*r^2 - 4*r^3/3;
   dwdr  = -4 + 8*r -4*r^2;
  
end
function [w,dwdr] = Spline5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 else
   w     = 1-10*r^3+15*r^4-6*r^5;
   dwdr  = -30*r^2 + 60*r^3-30*r^4;
  
end

function [w,dwdr] = power_function(arfa,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
else
    a2 = arfa*arfa;
   r2 = r*r;
    w     = exp(-r2/a2);
   dwdr  = (-2*r/a2)*exp(-r2/a2);
  
end

function [w,dwdr] = CSRBF2(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(6+36*r+82*r^2+72*r^3+30*r^4+5*r^5);
	dwdr  = 11*r*(r+2)*(5*r^3+15*r^2+18*r+4)*(r-1)^5;

end
function [w,dwdr] = CSRBF1(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^4*(4+16*r+12*r^2+3*r^3);
	dwdr  = -4*(1-r)^3*(4+16*r+12*r^2+3*r^3)+(1-r)^4*(16+24*r+9*r^2);

end
function [w,dwdr] = CSRBF3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
elseif(r==0)
    w =  1/3;
    dwdr  =  0.0;
else
	w     = 1/3+r^2-4/3*r^3+2*r^2*log(r);
	dwdr  = 4*r-4*r^2+4*r*log(r);
end
function [w,dwdr] = CSRBF4(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
elseif(r==0)
    w =  1/15;
    dwdr  =  0.0;
   
else
	w     = 1/15+19/6*r^2-16/3*r^3+3*r^4-16/15*r^5+1/6*r^6+2*r^2*log(r);
	dwdr  = 25/3*r-16*r^2+12*r^3-16/3*r^4+r^5+4*r*log(r);

end

function [w,dwdr] = CSRBF5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(35*r^2+18*r+3);
	dwdr  =-6*(1-r)^5*(35*r^2+18*r+3)+(1-r)^6*(70*r+18);
   
end

function [w,dwdr] = CSRBF6(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^8*(32*r^3+25*r^2+8*r+1);
	dwdr  =-8*(1-r)^7*(32*r^3+25*r^2+8*r+1)+(1-r)^8*(96*r^2+50*r+8);
   
end