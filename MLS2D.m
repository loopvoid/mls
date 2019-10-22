% two-DIMENSIONAL MLS APPROXIMATION
% by Jin Jia
%% 图像等步长采样
clc
clear all


I=imread('22result.jpg');
imshow(I);
[row,col,chn]=size(I);
% I=I(1:col,:,:);
% 设置节点坐标
step=10;%步长
xII=1: step : row;
yII=1: step : col;
[xI,yI] = meshgrid(yII,xII);
nnodes = size(xI,1)*size(yI,2);
% 设置评估点的坐标
[x,y] = meshgrid(1: 1 : col,1: 1: row);

npoints = size(x,1)*size(y,2);
%支持域范围
scale = 30;
% 确定每个节点的支持半径
dmI = scale *0.5* ones(1, nnodes);
tic
% 评估所有评估点x的MLS形状函数
[PHI, DPHIx, DPHIy] = MLS2DShape(6, nnodes, xI,yI, npoints, x,y, dmI, 'GAUSS', 3.0 ); 
toc
 
% 曲线拟合. y = peaks(x,y)
ZII  =I(xII,yII,:);    % 节点函数值 
% z  =x.*exp(-x.^2- y.^2);% 确切的解决方案
Zpoints=zeros(1,npoints);
xh=zeros(1,npoints);
yh=zeros(1,npoints);
II=I-I;
for i=1:npoints
% Zpoints(1,i)=z(i);
xh(1,i)=x(i);
yh(1,i)=y(i);
end
Znodes=zeros(1,nnodes);
for j=1:1
    ZI=ZII(:,:,j);
    for i=1:nnodes
        Znodes(1,i)=ZI(i);
    end                        %将二维数据转换为一维数据
    zh = PHI *Znodes';  % 逼近函数
    II(:,:,j)=reshape(zh,row,col);
end
ZI=double(ZI);
plot3( xI, yI, ZI,'k.','LineWidth',2);
hold on
surf(x,y,II(:,:,j));
toc
% imshow(II);
% III=imsubtract(I,II);
% sum(sum(sum(III)))/(row*col*3)
