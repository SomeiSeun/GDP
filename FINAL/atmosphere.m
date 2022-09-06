function [T,P,rho] = atmosphere(h)
% this funcction returns T in Kelvin and P in Pa for a given altitude
%   01/01/20 at Equator (0 long, 0 lat) at 00:00

L=length(h);

Data=load('atmosData.mat');


A=Data.atmosData(:,1);
B=Data.atmosData(:,2);
C=Data.atmosData(:,3);

T=zeros(L,1);
P=zeros(L,1);
for i=1:L
    [~,I]=min(abs(h(i)/1000-A));
    T(i,1)=C(I);
    if h > 500e3
        rho(i,1) = 0; 
        P(i,1) = 0;
    else
        rho(i,1)=1000*B(I);
        P(i,1)=rho(i)*287*T(i);
    end
end
end

