function [PHI, DPHIx, DPHIy] = MLS2DShape(m, nnodes, xI,yI, npoints, x,y, dmI, type, para)
% SHAPE FUNCTION OF 2D MLS APPROXIMATION
%
% SYNTAX: [PHI, DPHI, DDPHI] = MLS1DShape(m, nnodes, xI,yI, npoints, xi,yi, dmI, type, para)
%
% INPUT PARAMETERS
%    m - Total number of basis functions (1: Constant basis;  2: Linear basis;  3: Quadratic basis)
%    nnodes  - Total number of nodes used to construct MLS approximation
%    npoints - Total number of points whose MLS shape function to be evaluated
%    xI,yI(nnodes) - Coordinates of nodes used to construct MLS approximation
%    xi,yi(npoints) - Coordinates of points whose MLS shape function to be evaluated
%    dm(nnodes) - Radius of support of nodes
%    wtype - Type of weight function
%    para  - Weight function parameter
% ％m  - 基函数总数（1：常数基; 2：线性基; 3：二次基; 4）
% ％nnodes  - 用于构造MLS近似的节点总数
% ％npoints  - 要评估其MLS形状函数的点的总数
% ％xi（nnodes） - 用于构造MLS近似的节点坐标
% ％x（npoints） - 要评估其MLS形状函数的点的坐标
% ％dm（nnodes） - 节点支持的半径
% ％wtype  - 权重函数的类型
% ％para  - 权重函数参数
%
% OUTPUT PARAMETERS
%    PHI   - MLS Shpae function
%    DPHIx  - First order derivatives of MLS Shpae function to x
%    DPHIy - First order derivatives of MLS Shpae function to y
%
% INITIALIZE WEIGHT FUNCTION MATRICES
DmI=[];
wI   = zeros (nnodes, nnodes);  % Weight funciton
dwdxI  = zeros (1, nnodes);
dwdyI = zeros (1, nnodes);
xII = zeros(1,nnodes);
yII = zeros(1,nnodes);

% INITIALIZE SHAPE FUNCTION MATRICES
PHI   = zeros(npoints, nnodes);
DPHIx  = zeros(npoints, nnodes);
DPHIy = zeros(npoints, nnodes);

% LOOP OVER ALL EVALUATION POINTS TO CALCULATE VALUE OF SHAPE FUNCTION Fi(X)
for j = 1 : npoints
    DmI = dmI;
	% DETERMINE WEIGHT FUNCTIONS AND THEIR DERIVATIVES AT EVERY NODE
    [wII] =rectangleWeight(type, para, x(j),y(j),xI,yI,DmI,DmI);
    wI=diag(wII(:));
    xII=reshape(xI,size(xII));
    yII=reshape(yI,size(yII));
   % EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
   if (m == 1)  % Shepard function
      p = [ones(1, nnodes)]; 
      pxy   = [1];
      
      B    = p .* [wI];
     
   elseif (m == 3)
      B = [ones(1, nnodes); xII;yII]; 
      pxy   = [1; x(j);y(j)];
     
   elseif (m == 6)
      B = [ones(1, nnodes); xII;yII; xII.*xII;xII.*yII;yII.*yII]; 
      pxy   = [1; x(j); y(j);x(j)*x(j);x(j)*y(j);y(j)*y(j)];

   elseif (m == 10)
      B = [ones(1, nnodes); xII;yII; xII.*xII;xII.*yII;yII.*yII;xII.*xII.*xII;xII.*xII.*yII;xII.*yII.*yII;yII.*yII.*yII]; 
      pxy   = [1; x(j); y(j);x(j)*x(j);x(j)*y(j);y(j)*y(j);x(j)*x(j)*x(j);x(j)*x(j)*y(j);x(j)*y(j)*y(j);y(j)*y(j)*y(j)];
   else
      error('Invalid order of basis.');
   end
   
   % EVALUATE MATRICES A AND ITS DERIVATIVES
    B=B';
	A   = B'*wI*B;
   ARcond = rcond(A);
   
   
   while ARcond<=9.999999e-016  %判断条件数
       DmI=2*DmI;
       [wII] =rectangleWeight(type, para, x(j),y(j),xI,yI,DmI,DmI);
       wI=diag(wII(:));
       % EVALUATE MATRICES A AND ITS DERIVATIVES
        A   = B'*wI*B;
        ARcond = rcond(A);  
   end   %判断条件数   
      
   rxy  = pxy'/A;
   D=B'*wI;
   PHI(j,:) = rxy * D;   % shape function
    
end
