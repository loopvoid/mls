function [wi] = genGaussianWeightMat(dim,sigma)

% dim = 40;

h = fspecial('gaussian',dim,sigma);

in=0;
r=zeros(1,10);

for i = 1:dim
	for j=1:i
		in=in+1;
		r(1,in)=h(j,i);
	end
end

r_n=mapminmax(r,0,1);


out = zeros(dim,dim);
in=0;
for i=1:dim
	for j=1:i
		in=in+1;
		out(j,i)=r_n(1,in);
	end
end

out_t=out';

for i=1:dim
	for j=1:dim
		if i==j
			out_t(i,j)=0;
		end
	end
end

wi = out+out_t;
		