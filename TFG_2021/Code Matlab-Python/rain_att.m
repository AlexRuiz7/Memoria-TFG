clear all;
clc;
%variables
av=1.0085;
ah=1.1209;
kv=0.0001464;
kh=0.0001321;
theta=pi()*1/180*linspace(10,90,50);
tau=45*pi()/180;
R=2.89;
d_v=400;
d_hip=(400-.738)./sin(theta);
year_h=linspace(24, 24*365,365);
%calculo coeficients

for i=1:50
k(i)=(kh+kv+(kh-kv).*(cos(theta(i))).^2.*cos(2.*tau))./(2);
alpha(i)=(kh.*ah+kv.*av+(kh.*ah-kv.*av).*(cos(theta(i))).^2.*cos(2.*tau))./(2.*k(i));
end
Gamma=k.*(R.^alpha).*d_hip
plot(theta*180/pi(),Gamma)
title('Attenuation due to Rain with Elevation Angle')
xlabel('Elevation Angle (Degrees)');
ylabel('Attenuation (dB)');
%legend('Parte Real','Parte Imaginaria','Frec. con Impedancia ~ Real' ,'Frec. de diseño','Z=0');
