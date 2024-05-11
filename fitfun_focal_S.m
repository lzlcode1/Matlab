function [y] = fitfun_focal_S(x)

%%%%%%%%%%%%%%%%%     ！！！fitness function！！！   %%%%%%%%%%%%%%%%%%%%%%


freq = 15E6;%%%%Hz         ！！！Frequency setting！！！
d = 10E-3;%%%%m           ！！！Focal length setting！！！
L = 2E-2;%%%%m;           ！！！Lens radius！！！
a = 20E-6;%%%%m;          ！！！狭缝宽度！！！

c0 = 1500;      %        
rho0 = 1000;      %      
S0 = (a/2)^(3/2);
pixel_n = L/a;
pixel_n=round(pixel_n);
lambda=c0/freq;
k=2*pi/lambda;

m_v = x;
v = [m_v,fliplr(m_v)];
%v = zeros(1,pixel_n);
%v(1,100) = 1;
%v(1,500) = 1;

Rx = (0:a:a*(pixel_n-1))+a/2;
%R1 = repmat(Rx',1,pixel_n);
%R2 = repmat(Rx,pixel_n,1);
%R = abs(R1-R2);

%R_p = exp(-1i*k*R);
%R_f = R-eye(pixel_n)*(-2i*k*rho0*c0*S0);  
%func_G = R_p./sqrt(R_f);
%Zxr = -2i*k*rho0*c0*func_G*S0;
M_th = diag(v);
%Zxr = M_th*Zxr;
%Zxr = Zxr-diag(diag(Zxr));
Zxr = -rho0*c0*M_th+zeros(pixel_n);
Zxr = Zxr+eye(pixel_n)-M_th;

u0 = 1E0;
u_in = ones(pixel_n,1)*u0;

R_resis = diag(v);
intm = eye(pixel_n)-R_resis;
R_resis = R_resis*rho0*c0;
R_resis = R_resis+intm;

I = eye(pixel_n);
u_tran = (1/2)*(I-Zxr\R_resis)*u_in;

zoom = 10;
Rix = (0:a/zoom:L-a/zoom);
Ri1 = repmat(Rix',1,pixel_n);
Ri2 = repmat(Rx,pixel_n*zoom,1);
Rs = abs(Ri1-Ri2);

R_r = sqrt(Rs.^2+d^2); 
delta = a*sqrt(1-d./R_r)/2; 
func_G = exp(-1i*k*(R_r+delta))./sqrt(R_r);
Zst = 2i*rho0*c0*func_G*S0*k.*(sin(k*delta)./(k*delta+eps));
p_tran = Zst*u_tran;
pr = abs(p_tran);
pr = mapminmax(pr',0,1);

%plot(Rix,pr);
width = 100E-6;
c = width/a*zoom;
y1 = zeros(pixel_n*zoom,1);
y1(pixel_n/2*zoom-(c/2-1):pixel_n/2*zoom+1+(c/2-1)) = 50;
y = ones(1,pixel_n*zoom)*abs(y1-pr');

end