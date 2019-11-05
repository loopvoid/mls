% ?function least_square_fitting_Hyperbolic_paraboloid()

%% 程序说明%%

%b本程序是使用双曲抛物来进行

%验证最小二乘法，并且加入10%的噪音

%% figure 双曲抛物
x = -19:20;
y = -19:20;
[X,Y] =meshgrid(x,y);

Z = X.^2 - Y.^2;%准确解

figure
% mesh(X,Y,Z);
surf(X,Y,Z);
title('准确解')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% 加入10%噪声后进行拟合%%
A = zeros(6,6);%系数矩阵
b = zeros(1,6);
[row,col] =size(X);
        for i = 1:row
            for j = 1:col
                 A(1,1) = A(1,1) + 1;
                 A(1,2) = A(1,2) + X(i,j);
                 A(1,3) = A(1,3) + Y(i,j);
                 A(1,4) = A(1,4) + X(i,j)^2;
                 A(1,5) = A(1,5) +X(i,j)*Y(i,j);
                 A(1,6) = A(1,6) + Y(i,j)^2;

                 A(2,1) = A(2,1) + X(i,j);
                 A(2,2) = A(2,2) +X(i,j)*X(i,j)*X(i,j);
                 A(2,3) = A(2,3) +Y(i,j)*X(i,j);
                 A(2,4) = A(2,4) + X(i,j)^2*X(i,j);
                 A(2,5) = A(2,5) +X(i,j)*Y(i,j)*X(i,j);
                 A(2,6) = A(2,6) +Y(i,j)^2*X(i,j);

                 A(3,1) = A(3,1) + Y(i,j);
                 A(3,2) = A(3,2) + X(i,j)*Y(i,j);
                 A(3,3) = A(3,3) +Y(i,j)*Y(i,j);
                 A(3,4) = A(3,4) +X(i,j)^2*Y(i,j);
                 A(3,5) = A(3,5) +X(i,j)*Y(i,j)*Y(i,j);
                 A(3,6) = A(3,6) +Y(i,j)^2*Y(i,j);

                 A(4,1) = A(4,1) + X(i,j)^2;
                 A(4,2) = A(4,2) +X(i,j)*X(i,j)^2;
                 A(4,3) = A(4,3) +Y(i,j)*X(i,j)^2;
                 A(4,4) = A(4,4) +X(i,j)^2*X(i,j)^2;
                 A(4,5) = A(4,5) +X(i,j)*Y(i,j)*X(i,j)^2;
                 A(4,6) = A(4,6) +Y(i,j)^2*X(i,j)^2;

                 A(5,1) = A(5,1) +X(i,j)*Y(i,j);
                 A(5,2) = A(5,2) +X(i,j)*X(i,j)*Y(i,j);
                 A(5,3) = A(5,3) +Y(i,j)*X(i,j)*Y(i,j);
                 A(5,4) = A(5,4) + X(i,j)^2*X(i,j)*Y(i,j);
                 A(5,5) = A(5,5) +X(i,j)*Y(i,j)*X(i,j)*Y(i,j);
                 A(5,6) = A(5,6) +Y(i,j)^2*X(i,j)*Y(i,j);

                 A(6,1) = A(6,1) + Y(i,j)^2;
                 A(6,2) = A(6,2) +X(i,j)*Y(i,j)^2;
                 A(6,3) = A(6,3) +Y(i,j)*Y(i,j)^2;
                 A(6,4) = A(6,4) +X(i,j)^2*Y(i,j)^2;
                 A(6,5) = A(6,5) +X(i,j)*Y(i,j)*Y(i,j)^2;
                 A(6,6) = A(6,6) +Y(i,j)^2*Y(i,j)^2;
                 
                 b(1) = b(1) + Z(i,j);
                 b(2) = b(2) + Z(i,j)*X(i,j);
                 b(3) = b(3) + Z(i,j)*Y(i,j);
                 b(4) = b(4) + Z(i,j)*X(i,j)^2;
                 b(5) = b(5) +Z(i,j)*X(i,j)*Y(i,j);
                 b(6) = b(6) + Z(i,j)*Y(i,j)^2;
            end
        end

        a = b*inv(A);%求解系数
        x = X + 0.5*rand(row,col);%加入10%的噪音
        y = Y + 0.5*rand(row,col);%加入10%的噪音
        z = a(1) + a(2)*x + a(3)*y + a(4)*x.^2+ a(5)*x.*y + a(6)*y.^2;%拟合解

        figure
        mesh(X,Y,z)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title('加入10%噪音拟合结果')