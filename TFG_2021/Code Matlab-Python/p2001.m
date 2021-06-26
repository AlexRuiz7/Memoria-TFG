%Still a Work in Progress


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Header
%ITU-R P.2001-2,2015 RF Propagation Model
%A General Purpose Wide-Range Terrestrial Propagation Model in the Frequestion rang 30 MHz to 50 Ghz
%Distance: 3-1,000km
%Antenna Altitudes: Up to 8000m above sea level

function [Lb,debug_table] = p2001(freq,Tpol,rx_pos,tx_pos,Hrg,Htg,Tpc,Gr,Gt,h)
%cd('Z:\MATLAB\P2001')
%pause (0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Intro
%Section 2.2 Other Inputs:
%freq: frequency (GHz): 30 MHz - 50 GHz
%Tpol: Antenna Polarization: Vertical (1), Horizontal (0)
%rx_pos: Longitude, Latitude of Receiver (degrees)
%tx_pos: Longitude, Latitude of Trasnmitter (degrees)
%Hrg: Height of Electrical Center of Receiving Antenna above ground (meters)
%Htg: Height of Electrical Center of Transmitting Antenna above ground (meters)
%Tpc: Percentage of average year for which the prdicted basic transmission loss is not exceeded (%): Example 0.00001% to 99.99999%
%Gr: Gain (dBi) of receiving antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.
%Gt: Gain (dBi) of transmitting antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.
%h: Elevation data (evenly spaced) in meters.

%Constants: Section 2.3
lightspeed=physconst('lightspeed'); %Speed of Propagation: 299792458 m/s
Re=6371; %Average Earth Radius (km)
erland=22.0;  %Relative Permittivity for Land
ersea=80.0; %Relative Permittivity fo Sea
sigmaland=0.003; %Conducitivty for land (S/m)
sigmasea=5.0; %Conductivity for sea (S/m)


%Nick Reformat Data from .txt to .csv
DN_Median=csvread('nick_import_DN_Median.csv');
DN_SupSlope=csvread('nick_import_DN_SupSlope.csv');
DN_SubSlope=csvread('nick_import_DN_SubSlope.csv');
dndz_01=csvread('nick_import_dndz_01.csv'); 
TropoClim=csvread('nick_import_TropoClim.csv'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3
%Section 3.1 Limited Percentage Time
p_time=Tpc+0.00001*((50-Tpc)/50); %Equation 3.1.1
q_time=100-p_time; %Equation 3.1.2

%Great Circle Calculations: Attachment H or Matlab
dn=deg2km(distance(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2))); %Distance between tx and rx

FlagShort=dn<.1; %If the distance is less than 100m or 0.1km.

ellipsoid=[];
mid_npts=3;% This is based upon the number of points in the terrain model (h).
[t_mlat,t_mlon]=track2(rx_pos(1),rx_pos(2),tx_pos(1),tx_pos(2),ellipsoid,'degrees',mid_npts);
mid_point=[t_mlat(2),t_mlon(2)];

mid_npts5=5;% This is based upon the number of points in the terrain model (h).
[t2_mlat,t2_mlon]=track2(rx_pos(1),rx_pos(2),tx_pos(1),tx_pos(2),ellipsoid,'degrees',mid_npts5);

phi3qe=t2_mlon(2);
phi3qn=t2_mlat(2);

phi1qe=t2_mlon(4);
phi1qn=t2_mlat(4);

npts=length(h);% This is based upon the number of points in the terrain model (h).
[t_lat,t_lon]=track2(rx_pos(1),rx_pos(2),tx_pos(1),tx_pos(2),ellipsoid,'degrees',npts);
d=deg2km(distance(rx_pos(1),rx_pos(2),t_lat(1:length(t_lat)),t_lon(1:length(t_lat))));

mid_npts2=length(h);% This is based upon the number of points in the terrain model (h).
[tn_lat,tn_lon]=track2(rx_pos(1),rx_pos(2),tx_pos(1),tx_pos(2),ellipsoid,'degrees',mid_npts2);
for i=1:1:length(tn_lat)
    tropo_value(i)=get_txt_value(tn_lat(i),tn_lon(i),TropoClim);   %This was a guess as to the lat/lon input.
     %Gets the Climate of each point
end

%At this point, calculate dtm/dml on tropo_value map, only worrying about
%if it over land (1-6) or the sea (0) and using elevation data.
%Try to get access to the ITU Digitized World Map based on recommendation.

[dtm,dlm,dct,dcr,omega] = find_ITU_dist(rx_pos,tx_pos,h,dn);
h1=h(1);
hn=h(length(h));
profile_pts=length(h);

fsea=omega;
%Calculate Ground Height in meters above sea level at the mid-point
hmid=h(ceil(length(h)/2)); 
%hmid=h(446); %This was used for the validation where hmid=h(446)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.3 Antenna Altitude and Path Inclination
hts=h(1)+Htg; %Equation 3.3.1a 
hrs=h(length(h))+Hrg;  %Equation 3.3.1b

%Assign the higher and lower antenna height above sea level
hhi=max([hts,hrs]); %Equation 3.3.2a
hlo=min([hts,hrs]); %Equation 3.3.2b

%Calculate the Positive Value of Path Inclination:
ep=(hhi-hlo)/dn; %Equation 3.3.3
Sp=ep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.3
phitn=tx_pos(1);
phite=tx_pos(2);

phirn=rx_pos(1);
phire=rx_pos(2);

phime=mid_point(2); %Longitude Point at Midpoint
phimn=mid_point(1); %Latitude Point at Midpoint

Nd65m1=get_txt_value(phimn,phime,dndz_01); 

Sdn=get_txt_value(phimn,phime,DN_Median); %Sdn is interpolated from DN_Median.txt for the path mid-point
Nd1km50=-Sdn; %Equations 3.4.1.1 

SdeltaNsup=get_txt_value(phimn,phime,DN_SupSlope);
SdeltaNsub=get_txt_value(phimn,phime,DN_SubSlope);

if p_time<50
    Nd1kmp=Nd1km50+SdeltaNsup*log10(0.02*p_time); %Equation 3.4.1.2a
end
if p_time>=50
    Nd1kmp=Nd1km50-SdeltaNsub*log10(0.02*q_time); %Equation 3.4.1.2b
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.5 Effective Earth-radius geometry
ae=(157*Re)/(157+Nd1km50); %Equation 3.5.1 (km), Median Effective Earth Radius
Reff50=ae;
cp=(157+Nd1kmp)/(157*Re);   %Equation 3.5.2 (km-1)
%Effect Earth Curvature (cp is often positive, but it can be zero or negative.

%Effective Earth Radius exceeded for p% time limited but to become infinite
if cp>(10^-6)
    ap=1/cp; %Equation 3.5.3a
else
    ap=10^6; %Equation 3.5.3b
end
Reffp=ap;

%The path length expressed as the angle subtended by d km at the 
%center of a sphere of effective Earth radius:
thetae=dn/ae; %Original Equation 3.5.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.6 Wavelength
lambda=(10^-9*lightspeed)/freq; %Equation 3.6.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.7: Path Classification and terminal horizon parameters
%Find the highest elevation angle to an intermediate profile point, 
%relative to the horzontal at the transmitter.
clear temp_thetatim;
for i=2:1:(length(d)-1)
    temp_thetatim(i)=((h(i)-hts)/d(i))-((500*d(i))/ae); %Equation 3.7.1
end
[thetatim, ithetatim]=max(temp_thetatim(2:length(temp_thetatim))); %Equation 3.7.1
ithetatim=ithetatim+1;%To compensate for starting at 2

thetatr=((hrs-hts)/dn)-((500*dn)/ae); %Equation 3.7.2

%Case 1: Path is LOS 
FlagLos50=(thetatim<thetatr);
if thetatim<thetatr %The path is LOS
    
    for i=2:1:(length(d)-1)
        temp_vmax(i)=(h(i)+(((500*d(i))*(dn-d(i)))/ae)-((hts*(dn-d(i))+hrs*d(i))/dn))*sqrt((0.002*dn)/(lambda*d(i)*(dn-d(i)))); %Equation 3.7.3 
    end
    [vmax, ivmax]=max(temp_vmax(2:length(temp_vmax))); %Equation 3.7.3 
    ivmax=ivmax+1; %To compensate for starting at 2
    
    %Find the itermediate profile point with the highest diffraction parameter
    dim=d(ivmax);  %distance where imax is the profile index which gives vmax
    dlt=dim;  %Equation 3.7.4a (km)
    dlr=dn-dim; %Equation 3.7.4b (km)
    ilt=ivmax; %3.7.4c  
    ilr=ivmax; %3.7.4d
    thetat=thetatr; %Equation 3.7.5a (mrad)
    thetar=-thetatr-((1000*dn)/ae); %Equation 3.7.5b (mrad)
end

%Case 2: Path is NLOS
if thetatim>=thetatr %The path is NLOS
    dlt=d(ithetatim); %Equation 3.7.6a
    ilt=ithetatim; %Equation 3.7.6b
    thetat=thetatim; %Equation 3.7.7
    
    for i=2:1:(length(d)-1)
        temp_thetarim(i)=((h(i)-hrs)/(dn-d(i)))-((500*(dn-d(i)))/ae);%Equation 3.7.8
    end
    [thetarim, ithetarim]=max(temp_thetarim(2:length(temp_thetarim))); %Equation 3.7.8
    ithetarim=ithetarim+1; %To compensate for starting at 2
    dlr=dn-d(ithetarim);%Equation 3.7.9a
    ilr=ithetarim;%Equation 3.7.9b
    thetar=thetarim; %Equation 3.7.10
end

thetatpos=max([thetat, 0]); %Equation 3.7.11a
thetarpos=max([thetar, 0]); %Equation 3.7.11b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.8 Effective Heights and Path Roughness Parameter
v1=0;
for i=2:1:length(d)
   v1=v1+((d(i)-d(i-1))*(h(i)+h(i-1))); %Equation 3.8.1
end

v2=0;
for i=2:1:length(d)
   v2=v2+(d(i)-d(i-1))*(h(i)*(2*d(i)+d(i-1))+h(i-1)*(d(i)+2*d(i-1))); %Equation 3.8.2
end

hstip=(2*v1*dn-v2)/(dn*dn); %Equation 3.8.3a
hsrip=(v2-v1*dn)/(dn*dn); %Equation 3.8.3b

hstipa=min([hstip h(1)]); %Equation 3.8.4a
hsripa=min([hsrip h(length(h))]); %3.8.4.b

mses=(hsripa-hstipa)/dn; %Equation 3.8.5

htea=hts-hstipa; %Equation 3.8.6a
hrea=hrs-hsripa; %Equation 3.8.6b;

for i=ilt:1:ilr
    temp_hm(i)=h(i)-(hstipa+mses*d(i));%Equation 3.8.7
end
[hm, ihm]=max(temp_hm(ilt:1:ilr));%Equation 3.8.7

for i=2:1:(length(d)-1)
    capH(i)=h(i)-((hts*(dn-d(i))+hrs*d(i))/dn); %Equation 3.8.8d
end

hobs=max(capH(2:length(capH))); %Equation 3.8.8a

for i=2:1:(length(d)-1)
    temp_aobt(i)=capH(i)/d(i);%Equation 3.8.8b
end
aobt=max(temp_aobt(2:length(temp_aobt))); %Equation 3.8.8b

for i=2:1:(length(d)-1)
    temp_aobr(i)=capH(i)/(dn-d(i)); %Equation 3.8.8c
end
aobr=max(temp_aobr(2:length(temp_aobr))); %Equation 3.8.8c

if hobs<=0
    hst=hstip; %Equation 3.8.9a
    hsr=hsrip; %Equation 3.8.9b
else
    gt=aobt/(aobt+aobr); %Equation 3.8.9e
    gr=aobr/(aobt+aobr); %Equation 3.8.9f
    hst=hstip-hobs*gt; %Equation 3.8.9c
    hsr=hsrip-hobs*gr; %Equation 3.8.9d
end

if hst>h(1)
    hst=h(1);%Equation 3.8.10a
end

if hsr>h(length(h))
    hsr=h(length(h)); %Equation 3.8.10b
end

htep=hts-hst; %Equation 3.8.11a
hrep=hrs-hsr; %Equation 3.8.11b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.9 Tropospheric-scatter path segments
dtcv=(dn*tan(0.001*thetarpos+0.5*thetae)-0.001*(hts-hrs))/(tan(0.001*thetatpos+0.5*thetae)+tan(0.001*thetarpos+0.5*thetae)); %Equation 3.9.1a

if dtcv>dn %Limit dtcv to 0<=dtcv<=dn
   dtcv=dn; 
end
if dtcv<0 %Limit dtcv to 0<=dtcv<=dn
    dtcv=0;
end

drcv=dn-dtcv; %Equation 3.9.1b

%Calculate phicve, phicvn by setting dpnt=dtcv in H.3.1 (Great Circle Path)(No need to use Appendix H)
clear arclen;
clear az;
[arclen, az]=distance(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2));
[phicvn, phicve]=track1(tx_pos(1),tx_pos(2),az,km2deg(dtcv),ellipsoid,'degrees',1);

hcv=hts+1000*dtcv*tan(0.001*thetatpos)+((1000*dtcv*dtcv)/(2*ae)); %Equation 3.9.2

%Calculate phitcve, phitcvn by setting dpnt=0.5*dtcv in H.3.1 (Great Circle Path)(No need to use Appendix H)
clear arclen;
clear az;
[arclen, az]=distance(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2));
[phitcvn, phitcve]=track1(tx_pos(1),tx_pos(2),az,km2deg(dtcv*0.5),ellipsoid,'degrees',1);

%Calculate phitcve, phitcvn by setting dpnt=dn-0.5*drcv in H.3.1 (Great Circle Path)(No need to use Appendix H)
clear arclen;
clear az;
[arclen, az]=distance(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2));
[phircvn, phircve]=track1(tx_pos(1),tx_pos(2),az,km2deg(dn-0.5*drcv),ellipsoid,'degrees',1);

clear arclen;
clear az;
[arclen, az]=distance(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2));
[phicvn, phicve]=track1(tx_pos(1),tx_pos(2),az,km2deg(dtcv),ellipsoid,'degrees',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 3.10 Gaseous Absorption on Surface Paths
[Aosur,Awsur,Awrsur,gamma0,gammaw,gammawr]=f2_att(freq,phimn,phime,hmid,hts,hrs,dn);

Agsur=Aosur+Awsur; %Equation 3.10.1 The total gaseous attenuation under non-rain conditions

%Section 3.11 Free-space basic transmission loss
LbfsD=92.44+20*log10(freq)+20*log10(dn); %%Equation 3.11.1
Lbfs=LbfsD; %Equation 3.11.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 4: Sub-models
%There are 4 sub-models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 4.1: Sub-model 1: Normal propagation close to the surface of the Earth
%Use C.2 (Subroutine I)->B.2->B.4

%Perform the preliminary rain/wet-snow calculation in C.2 with the following inputs
phie=phime; %Equation 4.1.1a
phin=phimn; %Equation 4.1.1b
hrainlo=hlo; %Equation 4.1.1c
hrainhi=hhi; %Equation 4.1.1d
drain=dn;    %Equation 4.1.2

[Fwvr,Q0ra,Pm,Gm,hRtop,Pr6,dr,kmod,alphamod,a1,b1,c1] = c2_calc_mergeC5(phin,phie,q_time,hrainlo,hrainhi,drain,Tpol,freq,hlo,hhi); %Section C.2 Preliminary Calculation

[Q0ca]=b2_Q0ca(Nd65m1,thetatim,thetatr,dn,ep,hlo,freq,phimn,dlt,thetat,hts,ilt,h,ilr,d,thetar,hrs,dlr); %B.2 and B.3, Calculate Q0ca

[A1]=I_Aiter_submodel1(Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1,q_time,Q0ca);

%Start of Section A: Diffraction Loss

[Ldsph]=a2_loss(ap,htep,hrep,freq,dn,lambda,Tpol,omega); %Section A.2: Sphereical Earth Diffraction Loss

[Ldba,FlagLospa,Ldbka] = a4_Ldba(d,h,cp,hts,hrs,dn,lambda); %Section A.4: Bullington Diffraction Loss for Actual Profile

[Ldbs,FlagLosps,Ldbks]=a5_Ldbs(d,dn,ap,htep,hrep,lambda); %A.5: Bullington Diffraction Loss for a notional smooth profile

Ld=Ldba+max([Ldsph-Ldbs, 0]); %Equation A.1.1 

Lbm1=Lbfs+Ld+A1+Fwvr*(Awrsur-Awsur)+Agsur; %Equation 4.1.4 %OUTPUT of SUBMODEL 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4.1: Sub-model 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 4.2: Sub-model 2 Anomalous Propagation
[Lba,Aac,Aad,Aat] = d_Lba(phimn,dlt,dlr,thetat,thetar,freq,ae,dn,hm,htea,hrea,p_time,q_time,dtm,dlm,dct,dcr,omega); %Section D

Lbm2=Lba+Agsur; %Equation 4.2.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 4.3: Sub-model 3:Troposcatter propagation
[Lbs,thetas] = e_Lbs(phicvn,phicve,thetae,thetat,thetar,ae,dn,p_time,freq,Gt,Gr,Lbfs); %Use Attachment E to calculate Lbs

%Perform the preliminary rain/wet-snow calculation in C.2 for the TRANSMITTER to common-vloume path segment with the following inputs
phie=phitcve; %Equation 4.3.1a
phin=phitcvn; %Equation 4.3.1b
hrainlo=hts; %Equation 4.3.1c
hrainhi=hcv;    %Equation 4.3.1d
drain=dtcv;    %Equation 4.3.1e

%Save the value of Fwvr calculated in Section C.2 and call it Fwvrtx
[Fwvrtx,Q0ra,Pm,Gm,hRtop,Pr6,dr,kmod,alphamod,a1,b1,c1] = c2_calc_mergeC5(phin,phie,q_time,hrainlo,hrainhi,drain,Tpol,freq,hlo,hhi); %Section C.2 Preliminary Calculation

A2t=I_Aiter_submodel3(Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1,q_time);

%Perform the preliminary rain/wet-snow calculation in C.2 for the RECEIVER to common-vloume path segment with the following inputs
phie=phircve; %Equation 4.3.3a
phin=phircvn; %Equation 4.3.3b
hrainlo=hrs; %Equation 4.3.3c
hrainhi=hcv;    %Equation 4.3.3d
drain=drcv;    %Equation 4.3.3e


%Save the value of Fwvr calculated in Section C.2 and call it Fwvrrx
[Fwvrrx,Q0ra,Pm,Gm,hRtop,Pr6,dr,kmod,alphamod,a1,b1,c1] = c2_calc_mergeC5(phin,phie,q_time,hrainlo,hrainhi,drain,Tpol,freq,hlo,hhi); %Section C.2 Preliminary Calculation

A2r=I_Aiter_submodel3(Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1,q_time);

A2=(A2t*(1+0.018*dtcv)+A2r*(1+0.018*drcv))/(1+0.018*dn); %Equation 4.3.6

[Aos,Aws,Awrs,Aorcv,Aotcv,Awrcv,Awrrcv,Awrtcv,Awtcv,drcv]=f3_gas_absorption(h,phitn,phite,thetatpos,dtcv,phirn,phire,thetarpos,drcv,freq); %Use F.3 to Calculate Gaseous Attenuations Due to Oxygen

Ags=Aos+Aws; %Equation 4.3.7
Lbm3=Lbs+A2+0.5*(Fwvrtx+Fwvrrx)*(Awrs-Aws)+Ags; %Equation 4.3.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 4.4: Sporadic-E

[Lbe,foEs1,foEs2,Gamma1,Gamma2,LbEs1,LbEs2,Lp1r,Lp1t,Lp2r,Lp2t] = g_sporadic(p_time,phimn,phime,dn,freq,ae,phi1qe,phi1qn,phi3qe,phi3qn,thetat,thetar,dlt,dlr); %Use Attachment G to calculate Lbe
Lbm4=Lbe; %Equation 4.4.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start of Section 5
%Start of Section 5.1, Combining Sub-Models 1 and 2
clear Lm;
Lm=min([Lbm1, Lbm2]);
Lbm12=Lm-10*log10(10^(-0.1*(Lbm1-Lm))+10^(-0.1*(Lbm2-Lm))); %Equation 5.1.1
%End of Section 5.1

%Start of Section 5.2, Combining Sub Models 1+2, 3 and 4
clear Lm;
Lm=min([Lbm12, Lbm3, Lbm4]);
Lb=Lm-5*log10(10^(-0.2*(Lbm12-Lm))+10^(-0.2*(Lbm3-Lm))+10^(-0.2*(Lbm4-Lm))); %Equation 5.2.1
%End of Section 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 5
Blank=nan(1);

debug_file=table(Tpol,freq,Gr,Gt,Hrg,Htg,phire,phirn,phite,phitn,Tpc,Blank,FlagLos50,FlagLospa,FlagLosps,Blank,FlagShort,A1,A2,A2r,A2t,Aac,Aad,Aat,Ags,Agsur,Aorcv,Aos,Aosur,Aotcv,Awrcv,Awrrcv,Awrs,Awrsur,Awrtcv,Aws,Awsur,Awtcv,Blank,cp,dn,dcr,dct,Blank,dlm,dlr,dlt,drcv,dtcv,dtm,foEs1,foEs2,fsea,Fwvr,Fwvrrx,Fwvrtx,Gamma1,Gamma2,gamma0,gammaw,gammawr,h1,hcv,hhi,hlo,hm,hmid,hn,hrea,hrep,hrs,hsrip,hsripa,hstip,hstipa,htea,htep,hts,Lb,Lba,LbEs1,LbEs2,Lbfs,Lbm1,Lbm2,Lbm3,Lbm4,Lbs,Ld,Ldba,Ldbka,Ldbks,Ldbs,Ldsph,Lp1r,Lp1t,Lp2r,Lp2t,mses,profile_pts,Nd1km50,Nd1kmp,Nd65m1,ilr,ilt,Blank,Blank,Blank,Blank,phi1qe,phi1qn,phi3qe,phi3qn,phicve,phicvn,phime,phimn,phircve,phircvn,phitcve,phitcvn,Q0ca,Reff50,Reffp,Sp,thetae,thetar,thetarpos,thetas,thetat,thetatpos,p_time,q_time,Blank,lambda);                          
temp_cell=table2cell(debug_file);
temp2_cell=[debug_file.Properties.VariableNames;temp_cell];
debug_table=cell2table(temp2_cell');
end



%P.2001-2-2015: Attachement A.2: Spherical-Earch Diffraction Loss
function [Ldsph] = a2_loss(ap,htep,hrep,freq,dn,lambda,Tpol,omega)

    dlos=sqrt(2*ap)*(sqrt(0.001*htep)+sqrt(0.001*hrep)); %Equation A.2.1
    if dn>=dlos 
        adft=ap;
        %Calculate diffraction Loss using A.3 to find Ldft
        Ldft=a3_loss(freq,adft,Tpol,htep,hrep,omega,dn);  %Subrountine for A.3
        Ldsph=Ldft;
    else
        csph=(htep-hrep)/(htep+hrep); %Equation A.2.2d
        msph=(250*dn^2)/(ap*(htep+hrep)); %Equation A.2.2e
        
        bsph=2*sqrt((msph+1)/(3*msph))*cos(pi/3+1/3*acos((3*csph/2)*sqrt((3*msph)/((msph+1)^3)))); %Equation A.2.2c
        
        d1=dn/2*(1+bsph); %Equation A.2.2.a
        d2=dn-d1; %Equation A.2.2b
        
        hsph=((htep-500*(d1^2/ap))*d2+(hrep-500*(d2^2/ap))*d1)/dn; %Equation A.2.2
        
        hrep=17.456*sqrt((d1*d2*lambda)/dn); %Equation A.2.3
        
        if hsph>hrep
            Ldsph=0;
        else
            aem=500*(dn/(sqrt(htep)+sqrt(hrep)))^2; %Equation A.2.4
            adft=aem;
            %Use A.3 to calculate Ldft
            Ldft=a3_loss(freq,adft,Tpol,htep,hrep,omega,dn);  %Subrountine for A.3
            if Ldft<0
                Ldsph=0;
            else
                Ldsph=(1-(hsph/hrep))*Ldft; %Equation A.2.5
            end
        end
    end
end



%P.2001-2-2015: Attachement A.3: First-Term Spherical Earth Diffraction Loss using A.3 to find Ldft
function [Ldft] = a3_loss(freq,adft,Tpol,htep,hrep,omega,dn) 
    erland=22.0;  %Relative Permittivity for Land
    ersea=80.0; %Relative Permittivity fo Sea
    sigmaland=0.003; %Conducitivty for land (S/m)
    sigmasea=5.0; %Conductivity for sea (S/m)
   
    %%%%%LAND
    Ldftland=a3_subLdft(adft,freq,erland,sigmaland,Tpol,htep,hrep,dn); %Calculate Ldft using A.3.2 and A.3.8 to get Ldftland
    
    %%%%SEA
    Ldftsea=a3_subLdft(adft,freq,ersea,sigmasea,Tpol,htep,hrep,dn); %Calculate Ldft using A.3.2 and A.3.8 to get Ldftsea

    Ldft=omega*Ldftsea+(1-omega)*Ldftland; %Equation A.3.1
end


%P.2001-2-2015: Attachement A.3:Equations A.3.2-A.3.8 First-Term Spherical Earth Diffraction Loss using A.3 to find Ldft
function [subLdft] = a3_subLdft(adft,freq,er,sigma,Tpol,htep,hrep,dn)
    KH=0.036*(adft*freq)^(-1/3)*((er-1)^2+(18*sigma/freq)^2)^(-1/4); %Equation A.3.2a
    KV=KH*(er^2+(18*sigma/freq)^2)^(1/2); %Equation A.3.2b
    if Tpol==1
        Beta=(1+1.6*KV^2+0.67*KV^4)/(1+4.5*KV^2+1.53*KV^4); %Equation A.3.3
    else
        Beta=(1+1.6*KH^2+0.67*KH^4)/(1+4.5*KH^2+1.53*KH^4); %Equation A.3.3
    end
    X=21.88*Beta*((freq)/(adft^2))^(1/3)*dn; %Equation A.3.4
    Yt=0.9575*Beta*((freq^2)/adft)^(1/3)*htep; %Equation A.3.5a
    Yr=0.9575*Beta*((freq^2)/adft)^(1/3)*hrep; %Equation A.3.5b
    if X>=1.6
        FX=11+10*log10(X)-17.6*X; %Equation A.3.6
    else
        FX=-20*log10(X)-5.6488*X^(1.425); %Equation A.3.6
    end
    
    Bt=Beta*Yt; %Equation A.3.7a
    if Bt>2
        GYt=17.6*(Bt-1.1)^(0.5)-5*log10(Bt-1.1)-8; %Equation A.3.7
    else
        GYt=20*log10(Bt+0.1*Bt^3); %Equation A.3.7
    end
       
    Br=Beta*Yr; %Equation A.3.7a
    if Br>2
        GYr=17.6*(Br-1.1)^(0.5)-5*log10(Br-1.1)-8; %Equation A.3.7
    else
        GYr=20*log10(Br+0.1*Br^3); %Equation A.3.7
    end
    
    %Limit G(Y) such that G(Y)>= 2+20*log10(K)
    if Tpol==1 %Vertical
        if GYt<(2+20*log10(KV))
            GYt=2+20*log10(KV);
        end
        if GYr<(2+20*log10(KV))
            GYr=2+20*log10(KV);
        end
    else %Horizontal
        if GYt<(2+20*log10(KH))
            GYt=2+20*log10(KH);
        end
        if GYr<(2+20*log10(KH))
            GYr=2+20*log10(KH);
        end
    end
    subLdft=-FX-GYt-GYr; %Equation A.3.8 %Output
end


%P.2001-2-2015: Attachement A.4: Bullington Diffraction loss for actual profile
function [Ldba,FlagLospa,Ldbka] = a4_Ldba(d,h,cp,hts,hrs,dn,lambda)
    clear temp_Stim;
    for i=2:1:(length(d)-1)
        temp_Stim(i)=(h(i)+(500*cp*d(i)*(dn-d(i)))-hts)/d(i); %Equation A.4.1 
    end
    clear Stim;
    Stim=max(temp_Stim(2:length(temp_Stim))); %Equation A.4.1 

    Str=(hrs-hts)/dn; %Equation A.4.2

    FlagLospa=(Stim<Str);

    if Stim<Str %Case 1: LOS 
        clear temp_va;
        for i=2:1:(length(d)-1)
        	temp_va(i)=(h(i)+500*cp*d(i)*(dn-d(i))-((hts*(dn-d(i))+hrs*d(i))/dn))*sqrt((0.002*dn)/(lambda*d(i)*(dn-d(i)))); %Equation A.4.3  
        end
        vb=nan(1);
        db=nan(1);
        Srim=nan(1);
        va=max(temp_va(2:length(temp_va))); %Equation A.4.3          
        Ldbka=knife_edge_Jv(va); %Equation A.4.4  
    end
    
    if Stim>=Str %Case 2: NLOS 
        clear temp_Srim;
        for i=2:1:(length(d)-1)
            temp_Srim(i)=(h(i)+500*cp*d(i)*(dn-d(i))-hrs)/(dn-d(i)); %Equation A.4.5 {Correct}
        end
        clear Srim;
        Srim=max(temp_Srim(2:length(temp_Srim))); %Equation A.4.5 {Correct}
        clear db;    
        db=(hrs-hts+(Srim*dn))/(Stim+Srim); %Equation A.4.6 {Correct}
        clear vb;
        va=nan(1);
        vb=(hts+(Stim*db)-((hts*(dn-db)+(hrs*db))/dn))*sqrt((0.002*dn)/(lambda*db*(dn-db))); %Equation A.4.7 {Correct}
        clear Ldbka;
        Ldbka=knife_edge_Jv(vb); %Equation A.4.8
    end
    Ldba=Ldbka+(1-exp((-1*Ldbka)/6))*(10+(0.02*dn)); %Equation A.4.9 {Correct}
end

%P.2001-2-2015: Attachement A.5: Bullington Diffraction Loss for a notional smooth profile
function [Ldbs,FlagLosps,Ldbks] = a5_Ldbs(d,dn,ap,htep,hrep,lambda)
    clear temp_Stim;
    for i=2:1:(length(d)-1)
        temp_Stim(i)=(500*(dn-d(i))/ap)-(htep/d(i)); %Equation A.5.1
    end
    clear Stim;
    Stim=max(temp_Stim(2:length(temp_Stim))); %Equation A.5.1

    clear Str;
    Str=(hrep-htep)/dn; %Equation A.5.2 

    FlagLosps=(Stim<Str);
    if Stim<Str %Case 1: LOS 
        for i=2:1:(length(d)-1)
            vs1(i)=(500*d(i)*(dn-d(i)))/ap;
            vs2(i)=(htep*(dn-d(i))+(hrep*d(i)))/dn;
            vs3(i)=sqrt((0.002*dn)/(lambda*d(i)*(dn-d(i))));
        end
        
        for i=2:1:(length(d)-1)
           vs4(i)=(vs1(i)-vs2(i))*vs3(i);
        end
        vs=max(vs4(2:length(vs4))); %Equation A.5.3
        db=nan(1);
        vb=nan(1);
            
        Ldbks=knife_edge_Jv(vs); %Equation A.5.4
    end
    if Stim>=Str %Case 2: NLOS 
        
        clear temp_Srim;
        for i=2:1:(length(d)-1)
            temp_Srim(i)=((500*d(i))/ap)-(hrep/(dn-d(i))); %Equation A.5.5
        end
        clear Srim;
        Srim=max(temp_Srim(2:length(temp_Srim))); %Equation A.5.5
        
        vs=nan(1);
        clear db;
        db=(hrep-htep+Srim*dn)/(Stim+Srim); %Equation A.5.6
        
        clear vb;
        vb=((htep+Stim*db)-((htep*(dn-db)+hrep*db)/dn))*sqrt((0.002*dn)/(lambda*db*(dn-db))); %Equation A.5.7
        
        Ldbks=knife_edge_Jv(vb); %Equation A.5.8
    end
    
    Ldbs=Ldbks+(1-exp(-Ldbks/6))*(10+0.02*dn); %Equation A.5.9
end


%P.2001-2-2015: Attachement B.2: Characterize Multi-Path Activity
function [Q0ca,K,Q0cat,Q0car] = b2_Q0ca(Nd65m1,thetatim,thetatr,dn,ep,hlo,freq,phimn,dlt,thetat,hts,ilt,h,ilr,d,thetar,hrs,dlr)
    K=10^(-1*(4.6+0.0027*Nd65m1)); %Equation B.2.1
    %Case 1: Path is LOS 
    if thetatim<thetatr %The path is LOS
        dca=dn; %Equation B.2.2a
        eca=ep; %Equation B.2.2.b
        hca=hlo; %Equation B.2.2.c
        Q0car=nan(1);
        Q0cat=nan(1);
        Q0ca=b3_Q0ca(dca,eca,hca,freq,phimn,K); %USE B.3
    end

    %Case 2: Path is NLOS
    if thetatim>=thetatr %The path is NLOS
        dcat=dlt; %Equation B.2.3a
        ecat=abs(thetat); %Equation B.2.3b
        hcat=min([hts h(ilt)]); %Equation B.2.3c
    
        Q0cat=b3_Q0ca(dcat,ecat,hcat,freq,phimn,K); %USE B.3 to Calculate Q0cat

        %Calculate Q0car
        dcar=dlr; %Equation B.2.4a
        %dcar=46.348
        ecar=abs(thetar); %Equation B.2.4b
        hcar=min([hrs h(ilr)]); %Equation 4.2c
        
        Q0car=b3_Q0ca(dcar,ecar,hcar,freq,phimn,K); %USE B.3 to Calculate Q0car
        Q0ca=max([Q0cat Q0car]); %Equation B.2.5
    end
end


%P.2001-2-2015: Attachement B.3: Calculation of the Notional Zero-Fade Annual Percentage Time
function [Q0ca] = b3_Q0ca(dca,eca,hca,freq,phimn,K)
    qw=K*(dca^(3.1))*((1+eca)^(-1.29))*(freq^(0.8))*(10^(-0.00089*hca)); %Equation B.3.1
    if abs(phimn)<=45
        Cg=10.5-(5.6*log10(1.1+(abs(cosd(2*phimn))^(0.7))))-(2.7*log10(dca))+(1.7*log10(1+eca)); %Equation B.3.2a
    else
        Cg=10.5-(5.6*log10(1.1-(abs(cosd(2*phimn))^(0.7))))-(2.7*log10(dca))+(1.7*log10(1+eca)); %Equation B.3.2b
    end
    
    if Cg>10.8
        Cg=10.8;
    end
    Q0ca=10^(-0.1*Cg)*qw; %Equation B.3.3
end


%P.2001-2-2015: Attachement B.4: Percentage time a given clear-air fade level is exceeded on a surface path
function [Qcaf] = b4_Qcaf(Q0ca,A)
    if A>=0
        qt=3.576-1.955*log10(Q0ca); %Equation B.4.1.b
        qa=2+(1+0.3*10^(-0.05*A))*(10^(-0.016*A))*(qt+4.3*(10^(-0.05*A)+(A/800))); %Equation B.4.1a
        Qcaf=100*(1-exp(-10^(-0.05*qa*A)*log(2))); %B.4.1 %ln(2)
    end

    if A<0
        qs=-4.05-2.35*log10(Q0ca); %Equation B.4.2b
        qe=8+(1+0.3*10^(0.05*A))*(10^(0.035*A))*(qs+12*(10^(0.05*A)-(A/800))); %Equation B.4.2a
        Qcaf=100*exp(-10^(0.05*qe*A)*log(2)); %Equation B.4.2  %ln(2)
    end
end


%P.2001-2-2015: Attachement B.5: Percentage of time a given clear-air fade level is exceeded on a troposcatter path
function [Qcaftropo] = b5_Qcaftropo(A)
    if A<0
        Qcaftropo=100; %Equation B.5.1a
    else
        Qcaftropo=0; %Equation B.5.1b
    end
end

%P.2001-2-2015: Attachement A.2: Spherical-Earch Diffraction Loss
function [Fwvr,Q0ra,Pm,Gm,hRtop,Pr6,dr,kmod,alphamod,a1,b1,c1] = c2_calc_mergeC5(phin,phie,q_time,hrainlo,hrainhi,drain,Tpol,freq,hlo,hhi)

    %Get Pr6 from 'Esarain_Pr6_v5.txt' using phin and phie
    %Get MT from 'Esarain_MT_v5.txt' using phin and phie
    %Get Brain from 'Esarain_Beta_v5.txt' using phin and phie

    Esarain_Pr6_v5=csvread('nick_import_Esarain_Pr6_v5.csv');
    Esarain_Mt_v5=csvread('nick_import_Esarain_Mt_v5.csv');
    Esarain_Beta_v5=csvread('nick_import_Esarain_Beta_v5.csv');
    h0text=csvread('nick_import_h0.csv');

    Pr6=get_txt_value(phin,phie,Esarain_Pr6_v5);
    MT=get_txt_value(phin,phie,Esarain_Mt_v5);
    Brain=get_txt_value(phin,phie,Esarain_Beta_v5);
    h0=get_txt_value(phin,phie,h0text);

    hR=360+1000*h0; %Equation C.2.1
    hRtop=hR+2400; %Equation C.2.2

    if Pr6==0 || hrainlo>=hRtop %No Rain,
        Q0ra=0; %Main outputs
        Fwvr=0; %Main outputs

        Pm=0; %Place Holders when there is no rain
        Gm=0; %Place Holders when there is no rain
        dr=0; %Place Holders when there is no rain
        kmod=0; %Place Holders when there is no rain
        alphamod=0; %Place Holders when there is no rain
        a1=0; %Place Holders when there is no rain
        b1=0; %Place Holders when there is no rain
        c1=0; %Place Holders when there is no rain
    else
        Mc=Brain*MT; %Equation C.2.3a {Correct}
        MS=(1-Brain)*MT; %Equation C.2.3b
        Q0ra=Pr6*(1-exp((-0.0079*MS)/Pr6)); %Equation C.2.4
        a1=1.09; %C.2.5a
        b1=(Mc+MS)/(21797*Q0ra); %Equation C.2.5b
        c1=26.02*b1; %Equation C.2.5c
        Qtran=Q0ra*exp((a1*(2*b1-c1))/(c1*c1)); %Equation C.2.6

        if Tpol==1 %Vertical
            tau=90;
        else
            tau=0;
        end

        erain=(0.001*(hrainhi-hrainlo))/drain; %Equation C.2.7

        %Use ITU-R P.838 to calculate rain regression coefficients for k and alpha
        if freq<1
            [k1ghz, alpha1ghz]=rain_coefficients_838(1,deg2rad(tau),erain);
            k=freq*k1ghz; %Equation C.2.8a %k1ghz are calculated for 1GHz
            alpha=alpha1ghz; %Equation C.2.8b %alpha1ghz are calculated for 1GHz
        else
            [k, alpha]=rain_coefficients_838(freq,deg2rad(tau),erain);
        end

        dr=min([drain 300]); %Equation C.2.9a
        drmin=max([dr 1]); %Equation C.2.9b

        kmod1=1.763^(alpha)*k;
        kmod2=(0.6546*exp(-0.009516*drmin)+0.3499*exp(-0.001182*drmin));
        kmod=kmod1*kmod2; %Equation C.2.10a

        alphamod1=(0.753+(0.197/drmin))*alpha;
        alphamod2=0.1572*exp(-0.02268*drmin);
        alphamod3=0.1594*exp(-0.0003617*drmin);

        alphamod=alphamod1+alphamod2-alphamod3; %Equation C.2.10b

        Pm=zeros(1,49); %Proability of Particular case, initialize all Pm to zero
        Gm=zeros(1,49);
        Gm(1)=1;
        Hn=-2400:100:2400;
        m=1;  %marker

        %For each line of Table C.2.1, for n from 1 to 49, do the following:
        for n=1:1:length(Pm)
            %a)
            clear hT; %Equation C2. a)
            hT=hR+Hn(n); %Hn is the corresponding relative height entry in Table C.2.1
            %b)
            if hrainlo>=hT %Equation C2. b)
                %repaet from a) for the next value of n
                %Do nothing
            else
                %do c)
                if hrainhi>(hT-1200) %Equation C2. c)
                    %c.i)Use the method C.5 to set Gm to the path-averaged multiplier for this path geometry relative ot the melting layer
                    clear slo;
                    clear shi;
                    slo=1+floor((hT-hlo)/100); %Equation C.5.1a
                    shi=1+floor((hT-hhi)/100); %Equation C.5.1b
                    if slo<1
                        G=0;
                        %End of G calculation
                    elseif shi>12
                        G=1;
                        %End of G claculation
                    elseif slo==shi
                        clear deltah;
                        deltah=0.5*(hlo+hhi)-hT; %Equation C.5.2
                        GAMMAslice=c4_Gamma(deltah);  %Equation C.5.2
                        G=GAMMAslice;
                        %End of claculation
                    else
                        G=0; %initialize, Equation C.5.3
                        s(1)=max([shi,1]); %Equation C.5.4a %sfirst
                        s(2)=min([slo,12]); %Equation C.5.4b %slast

                        %for j=1:1:length(s)
                        for j=s(1):1:s(2)
                            %Condition 1
                            if shi<j && j<slo
                                deltah=100*(0.5-j); %Equation C.5.5a
                                Q=100/(hhi-hlo); %Equation C.5.5b
                            elseif j==slo %Condition 2
                                deltah=0.5*((hlo-hT-(100*(j-1)))); %Equation C.5.6a
                                Q=(hT-(100*(j-1))-hlo)/(hhi-hlo); %Equation C.5.6b
                            elseif j==shi %Condition 3
                                deltah=0.5*(hhi-hT-(100*j)); %Equation C.5.7a
                                Q=(hhi-(hT-100*j))/(hhi-hlo); %Equation C.5.7b
                            end

                            GAMMAslice=c4_Gamma(deltah); %Equation C.5.8;
                            G=G+(Q*GAMMAslice); %Equation C.5.9
                        end
                    end

                    if slo>12
                        clear Q;
                        Q=(hT-1200-hlo)/(hhi-hlo); %Equation C.5.10
                        G=G+Q; %Equation C.5.11
                    end
                    Gm(m)=G;
                    %c.ii) set Pm= to the third column (Cap pi sub n) in table C.2.1
                    Pm(m)=prob_C21(n);
                    %c.iii)
                    if n<49
                        m=m+1; %add 1 to array index m
                    end
                    %c.iv) repeat from a) for the next value of n

                else
                    %Continue to d)
                    %d) Accumulate (Cap pi sub n) from Table C.2.1 into Pm,
                    %   set Gm=1 and repeat a) for the next value of n.
                    if Pm(m)==0
                        Pm(m)=prob_C21(n);
                    else
                        clear Pm_temp;
                        Pm_temp=Pm(m);
                        Pm(m)=Pm_temp+prob_C21(n);
                    end
                    Gm(m)=1;
                end
            end
        end

        %Set the number of values in arrays Gm and Pm according to M=m Equation C.2.12
        clear temp_Pm
        temp_Pm=Pm(1:m);
        clear Pm;
        Pm=temp_Pm;

        clear temp_Gm
        temp_Gm=Gm(1:m);
        clear Gm;
        Gm=temp_Gm;
        clear Rwvr;
        Rwvr=6*((log10(Q0ra/q_time))/(log10(Q0ra/Qtran)))-3; %Equation C.2.13a
        clear temp_sum;
        for i=1:1:length(Gm)
            temp_sum(i)=Gm(i)*Pm(i);%Equation C.2.13pre
        end
        clear prod_GmPm;
        prod_GmPm=sum(temp_sum); %Equation C.2.13pre
        Fwvr=0.5*(1+tanh(Rwvr))*prod_GmPm; %Equation C.2.13 {Correct}
    end
end

%VALIDATED 
%ITU-R P.838-3 Rain Attenuation
function [k alpha]=rain_coefficients_838(freq,tau,erain)
    %Use ITU-R P.838 to calculate rain regression coefficients for k and alpha
    
    %Tables from ITU-R 838.3 (Tables 1-4)
    kH_aj=[-5.33980, -0.35351, -0.23789, -0.94158];
    kH_bj=[-0.10008, 1.26970, 0.86036, 0.64552];
    kH_cj=[1.13098,0.45400, 0.15354, 0.16817];
    kH_mk=-0.18961;
    kH_ck=0.71147;
    
    kV_aj=[-3.80595,-3.44965,-0.39902,0.50167];
    kV_bj=[0.56934,-0.22911,0.73042,1.07319];
    kV_cj=[0.81061,0.51059,0.11899,0.27195];
    kV_mk=-0.16398;
    kV_ck=0.63297;
    
    aH_aj=[-0.14318, 0.29591, 0.32177, -5.37610, 16.1721];
    aH_bj=[1.82442, 0.77564, 0.63773, -0.96230, -3.29980];
    aH_cj=[-0.55187, 0.19822, 0.13164, 1.47828, 3.43990];
    aH_ma=0.67849;
    aH_ca=-1.95537;
    
    aV_aj=[-0.07771, 0.56727, -0.20238, -48.2991, 48.5833];
    aV_bj=[2.33840, 0.95545, 1.14520, 0.791669, 0.791459];
    aV_cj=[-0.76284, 0.54039, 0.26809, 0.116226, 0.116479];
    aV_ma=-0.053739;
    aV_ca=0.83433;
    
    clear temp_kH;
    for j=1:1:4 
        temp_kH(j)=kH_aj(j)*exp(-1*((log10(freq)-kH_bj(j))/kH_cj(j))^2);  %ITU-R 838.3 Equation 2
    end
    sum_kH=sum(temp_kH)+kH_mk*log10(freq)+kH_ck; %ITU-R 838.3 Equation 2
    kH=10^sum_kH; %ITU-R 838.3 Equation 2
    
    clear temp_alphaH;
    for j=1:1:5 
        temp_alphaH(j)=aH_aj(j)*exp(-1*((log10(freq)-aH_bj(j))/aH_cj(j))^2);  %ITU-R 838.3 Equation 3
    end
    alphaH=sum(temp_alphaH)+aH_ma*log10(freq)+aH_ca; %ITU-R 838.3 Equation 3
    
    clear temp_kV;
    for j=1:1:4 
        temp_kV(j)=kV_aj(j)*exp(-1*((log10(freq)-kV_bj(j))/kV_cj(j))^2);  %ITU-R 838.3 Equation 2
    end
    sum_kV=sum(temp_kV)+kV_mk*log10(freq)+kV_ck; %ITU-R 838.3 Equation 2
    kV=10^sum_kV; %ITU-R 838.3 Equation 2
    
    clear temp_alphaV;
    for j=1:1:5 
        temp_alphaV(j)=aV_aj(j)*exp(-1*((log10(freq)-aV_bj(j))/aV_cj(j))^2);  %ITU-R 838.3 Equation 3
    end
    alphaV=sum(temp_alphaV)+aV_ma*log10(freq)+aV_ca; %ITU-R 838.3 Equation 3
          
     k=(kH+kV+(kH-kV)*cos(erain)^2*cos(2*tau))/2; %ITU-R 838.3 Equation 4
     alpha=(kH*alphaH+kV*alphaV+(kH*alphaH-kV*alphaV)*cos(erain)^2*cos(2*tau))/(2*k); %ITU-R 838.3 Equation 5
end


%P.2001-2-2015: Find Value of Table C.2.1
function [value] = prob_C21(index) 
probability_C21=[0.000555,
    0.000802,
    0.001139,
    0.001594,
    0.002196,
    0.002978,
    0.003976,
    0.005227,
    0.006764,
    0.008617,
    0.010808,
    0.013346,
    0.016225,
    0.019419,
    0.022881,
    0.026542,
    0.030312,
    0.034081,
    0.037724,
    0.041110,
    0.044104,
    0.046583,
    0.048439,
    0.049589,
    0.049978,
    0.049589,
    0.048439,
    0.046583,
    0.044104,
    0.041110,
    0.037724,
    0.034081,
    0.030312,
    0.026542,
    0.022881,
    0.019419,
    0.016225,
    0.013346,
    0.010808,
    0.008617,
    0.006764,
    0.005227,
    0.003976,
    0.002978,
    0.002196,
    0.001594,
    0.001139,
    0.000802,
    0.000555];
value=probability_C21(index);
end


%P.2001-2-2015: Attachement C.3: Percentage Time a Given Precipitation Fade Level is Exceeded
function [Qrain] = c3_Qrain(A,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1)
    if A<0 
        Qrain=100; %C.3.1a  
    end
    if A>=0 %This is dependent on whether the path is classifed as "rain" or "non rain"  
        if Pr6==0 || hrainlo>=hRtop %No Rain, %Pulled from C.2
            Qrain=0; %Equation C.3.1b, "Non-rain"
        else %"Rain", Where these Pm and Gm are defined in C.2.11-C.2.13a
            drlim=max([dr, 0.001]);  %Equation C.3.1e [km]
            clear temp_Qrain;
            clear Rm;
            for m=1:1:length(Gm)
                Rm(m)=(A/(Gm(m)*drlim*kmod))^(1/alphamod); %Equation C.3.1d
                temp_Qrain(m)=(Pm(m)*exp(-1*(a1*Rm(m)*(b1*Rm(m)+1))/(c1*Rm(m)+1))); %Equation C.3.1c
            end
            Qrain=100*sum(temp_Qrain);
        end
    end
end


%P.2001-2-2015: C.4: Melting Layer Model, returning the attenuation multiplier
function [capGamma] = c4_Gamma(deltah) 
    if 0<deltah
        capGamma=0;%Equation C.4.1
    elseif -1200<=deltah && deltah<=0
        gamma1=4*(1-exp(deltah/70))^2;
        gamma2=(1-exp(-((deltah/600)^2)))^2;
        gamma3=(4*((1-exp(deltah/70))^2)-1);
        capGamma=gamma1/(1+(gamma2*gamma3));
    elseif deltah<-1200
        capGamma=1;%Equation C.4.1
    end
end

%P.2001-2-2015: Attachment D: Anomalous/Layer-Reflection Model
function [Lba,Aac,Aad,Aat] = d_Lba(phimn,dlt,dlr,thetat,thetar,freq,ae,dn,hm,htea,hrea,p_time,q_time,dtm,dlm,dct,dcr,omega)
    %Start of D.2: Point Incidence of ducting
    tao=round(1-exp(-0.000412*(dlm^2.41)),6); %Equation D.2.1
    
    u1=(10^((-1*dtm)/(16-6.6*tao))+10^(-1*(2.48+1.77*tao)))^0.2; %Equation D.2.2
    
    if u1>1 %Limit the value of u1<=1
        u1=1;
    end
    
    if abs(phimn)<=70 
        u4=10^((-0.935+0.0176*abs(phimn))*log10(u1)); %Equation D.2.3
    else
        u4=10^(0.3*log10(u1)); %Equation D.2.3
    end
    
    if abs(phimn)<=70
        Beta0=(10^(-0.015*abs(phimn)+1.67))*u1*u4; %Equation D.2.4
    else
        Beta0=4.17*u1*u4; %Equation D.2.4
    end
    
    %End of D.2
    %Start of D.3: Site-Shielding losses with respect to anomalous propagation mechanisms
    gtr=0.1*dlt; %Equation D.3.1a
    grr=0.1*dlr; %Equation D.3.1b
    
    thetast=thetat-gtr; %Equation D.3.2a
    thetasr=thetar-grr; %Equation D.3.2b
    
    if thetast>0
        Ast=20*log10(1+0.361*thetast*(freq*dlt)^(1/2))+0.264*thetast*freq^(1/3); %Equation D.3.3a
    else
        Ast=0; %Equation D.3.3b
    end

    if thetasr>0
        Asr=20*log10(1+0.361*thetasr*(freq*dlr)^(1/2))+0.264*thetasr*freq^(1/3); %Equation D.3.4a
    else
        Asr=0;  %Equation D.3.4b
    end
    %End of D.3
    
    %Start of D.4: Over-sea surface duct coupling corrections
    if omega>=0.75 && dct<=dlt && dct<=5 %5km
        Act=-3*exp(-0.25*(dct^2))*(1+tanh(0.07*(50-hts))); %Equation D.4.2a
    else
        Act=0; %Equation D.4.2b
    end
    
    if omega>=0.75 && dcr<=dlr && dcr<=5 %5km
        Acr=-3*exp(-0.25*dcr^2)*(1+tanh(0.07*(50-hrs))); %Equation D.4.3a
    else
        Acr=0;  %Equation D.4.3b
    end
    
    %End of D.4
    
    %Start of D.5
    if freq<0.5
        Alf=(45.375-137*freq+92.5*freq^2)*omega; %Equation D.5.2a
    else
        Alf=0;  %Equation D.5.2b
    end
    
    Aac=102.45+20*log10(freq*(dlt+dlr))+Alf+Ast+Asr+Act+Acr; %Equation D.5.1
    
    %End of D.5
    
    %Start of D.6
    
    gammad=5*(10^-5)*ae*(freq^(1/3)); %Equation D.6.1
    
    thetaat=min([thetat, gtr]); %Equation D.6.2a
    thetaar=min([thetar, grr]); %Equation D.6.2b
    
    thetaa=((1000*dn)/ae)+thetaat+thetaar; %Equation D.6.3
    
    Aad=gammad*thetaa; %Equation D.6.4 (Correct EQ)
    
    %End of D.6
    
    %Start of D.7
    
    dar=min([dn-dlt-dlr, 40]); %Equation D.7.1
    
    if hm>10
        u3=exp(-4.6*(10^-5)*(hm-10)*(43+(6*dar))); %Equation D.7.2a
    else
        u3=1; %Equation D.7.2b
    end
    
    alpha=-0.6-3.5*(10^-9)*(dn^3.1)*tao; %Equation D.7.3
    
    if alpha<(-3.4)
        alpha=-3.4;
    end
    
    u2=((500*dn^2)/(ae*(sqrt(htea)+sqrt(hrea))^2))^alpha; %Equation D.7.4
    
    if u2>1
        u2=1;
    end
    
    Betaduct=Beta0*u2*u3; %D.7.5
    
    CapGamma=(1.076*exp((-10^-6)*(dn^1.13)*(9.51-4.8*log10(Betaduct)+0.198*(log10(Betaduct))^2)))/((2.0058-log10(Betaduct))^1.012); %Equation D.7.6
    
    Aat=-12+((1.2+0.0037*dn)*log10(p_time/Betaduct))+(12*(p_time/Betaduct)^CapGamma)+(50/q_time); %Equation D.7.7
    %End of D.7
    
    Lba=Aac+Aad+Aat; %Equation D.8.1
end

%P.2001-2-2015: Attachment E: Troposcatter
function [Lbs,thetas] = e_Lbs(phicvn,phicve,thetae,thetat,thetar,ae,dn,p_time,freq,Gt,Gr,Lbfs)
    TropoClim=csvread('nick_import_TropoClim.csv'); 
    mp_trop_cli=get_txt_value(phicvn,phicve,TropoClim);
    
    M=[129.60, 119.73, 109.30, 128.50, 119.73, 123.2, 116]; %From Table E.1
    gamma=[0.33, 0.27, 0.32, 0.27, 0.27, 0.27, 0.27]; %From Table E.1

    thetas=1000*thetae+thetat+thetar; %Equation E.1
    
    htrop=0.125*10^(-6)*thetas^2*ae; %Equation E.4
    He=0.25*10^(-3)*dn; %Equation E.3
   
    if mp_trop_cli==0
        LN=20*log10(5+gamma(7)*He)+4.34*gamma(7)*htrop; %Equation E.2
    else
        LN=20*log10(5+gamma(mp_trop_cli)*He)+4.34*gamma(mp_trop_cli)*htrop; %Equation E.2
    end
    
    ds=0.001*thetas*ae; %Equation E.5

    if mp_trop_cli==1
        if ds<100
            Y90=-8.2;%Equation E.8
        end
        if ds<1000 && ds>=100
            Y90=1.006*10^(-8)*(ds^3)-2.569*10^(-5)*(ds^2)+0.2242*ds-10.2; %Equation E.8
        end
        if ds>=1000
            Y90=-3.4; %Equation E.8
        end
    end
    if any(mp_trop_cli==[2,5,6])
        Y90=-2.2-(8.1-0.23*min([freq, 4]))*exp(-0.137*htrop); %Equation E.6
    end
    if mp_trop_cli==3
        if ds<100
            Y90=-10.845; %Equation E.9
        end
        if ds<465 && ds>=100
            Y90=-4.5*10^(-7)*(ds^3)+4.45*10^(-4)*(ds^2)-0.122*ds-2.645; %Equation E.9
        end
        if ds>=465
            Y90=-8.4; %Equation E.9
        end
    end
    if mp_trop_cli==4
        if ds<100
            Y90=-11.5;%Equation E.10
        end
        if ds<550 && ds>=100
            Y90=-8.519*10^(-8)*(ds^3)+7.444*10^(-5)*(ds^2)-4.18*10^(-4)*ds-12.1; %Equation E.10
        end
        if ds>=550
            Y90=-4.0; %Equation E.10
        end    
    end
    if mp_trop_cli==0 %'Sea'
        Y90=-9.5-3*exp(-0.137*htrop); %Equation E.7
    end

    if p_time>=50
        C=1.26*(-1*log10((100-p_time)/50))^(0.63); %Euqation E.11a
    else
        C=-1.26*(-1*log10(p_time/50))^0.63; %Equation E.11b
    end
    Yp=C*Y90; %Equation E.12

    if thetas<10^(-6) %Limit the value of thetas to thetas>=10^-6
        thetas=10^(-6);
    end
   
    Ldist=max([10*log10(dn)+30*log10(thetas)+LN, 20*log10(dn)+0.573*thetas+20]); %Equation E.13
    Lfreq=25*log10(freq)-2.5*(log10(0.5*freq))^2; %Equation E.14
    Lcoup=0.07*exp(0.055*(Gt+Gr)); %Equation E.15
    
    if mp_trop_cli==0
        Lbs=M(7)+Lfreq+Ldist+Lcoup-Yp; %Equation E.16
    else
        Lbs=M(mp_trop_cli)+Lfreq+Ldist+Lcoup-Yp; %Equation E.16
    end

    %To Avoid under-estimating troposcatter loss for short paths, limits Lbs
    if Lbs<Lbfs
       Lbs=Lbfs; %Equation E.17 
    end
end

%P.2001-2-2015: Attachment F2/F5/F6: Gaseous Absoprtion for Surface Path
%Not valid for frequencies greater than 54 GHz, but P2001 is limited to 50 GHz
%See ITU-R P.676 for a more general expression
function [Aosur,Awsur,Awrsur,gamma0,gammaw,gammawr] = f2_att(freq,phimn,phime,hmid,hts,hrs,dn)
    hsur=hmid; %The terrain heigh at the middle of the path
    [gamma0,gammaw,gammawr] = f6_att(freq,phimn,phime,hsur);
    hrho=0.5*(hts+hrs); %Equation F.2.1

    Aosur=gamma0*dn*exp(-hrho/5000); %Equation F.2.2a
    Awsur=gammaw*dn*exp(-hrho/2000); %Equation F.2.2b
    Awrsur=gammawr*dn*exp(-hrho/2000); %Equation F.2.2c
end

%P.2001-2-2015: Attachment F.3 Gaseous Absorption for a Troposcatter Path
function [Aos,Aws,Awrs,Aorcv,Aotcv,Awrcv,Awrrcv,Awrtcv,Awtcv,drcv] = f3_gas_absorption(h,phitn,phite,thetatpos,dtcv,phirn,phire,thetarpos,drcv,freq)
    hsur=h(1);
    thetaelev=thetatpos;
    dcv=dtcv;
    [Aotcv,Awtcv,Awrtcv] = f4_path(hsur,thetaelev,dcv,phitn,phite,freq); %Use Equation F.3.1a-c

    clear hsur;
    hsur=h(length(h));
    clear thetaelve;
    thetaelev=thetarpos;
    clear dcv;
    dcv=drcv;

    [Aorcv,Awrcv,Awrrcv] = f4_path(hsur,thetaelev,dcv,phirn,phire,freq); %Use Equation F.3.2a-c

    Aos=Aotcv+Aorcv; %Equation F.3.3a
    Aws=Awtcv+Awrcv; %Equation F.3.3b
    Awrs=Awrtcv+Awrrcv; %Equation F.3.3c
end


%P.2001-2-2015: Attachment F.4 Gaseous Absorption for a Terminal/Common-Volume Troposcatter Path
function [Ao,Aw,Awr] = f4_path(hsur,thetaelev,dcv,lat,lon,freq)
    [gamma0,gammaw,gammawr] = f6_att(freq,lat,lon,hsur); %Use F.6.2
    
    do=5/(0.65*sin(0.001*thetaelev)+0.35*sqrt((sin(0.001*thetaelev).^2)+0.00304)); %Equation F.4.1a
    dw=2/(0.65*sin(0.001*thetaelev)+0.35*sqrt((sin(0.001*thetaelev).^2)+0.00122)); %Equation F.4.1b

    deo=do*(1-exp(-1*dcv/do))*exp(-1*hsur/5000); %Equation F.4.2a
    dew=dw*(1-exp(-1*dcv/dw))*exp(-1*hsur/2000); %Equation F.4.2.b
    
    Ao=gamma0*deo; %Equation F.4.3a
    Aw=gammaw*dew; %Equation F.4.3b
    Awr=gammawr*dew; %Equation F.4.3c
end


%P.2001-2-2015: Attachment F6(F5): Gaseous Absoprtion for Surface Path
%Not valid for frequencies greater than 54 GHz. P2001 is limited to 50 GHz
%See ITU-R P.676 for a more general expression
function [gamma0,gammaw,gammawr] = f6_att(freq,lat,lon,hsur)
    %F.6 Specific Sea Level Attenuations
    gamma0=((7.2/(freq^2+0.34))+(0.62/((54-freq)^1.16+0.83)))*freq^2*10^-3; %Equation F.6.1

    Surfwv_50_fixed=csvread('nick_import_Surfwv_50_fixed.csv');
    psur=get_txt_value(lat,lon,Surfwv_50_fixed);

    psea=psur*exp(hsur/2000); %Equation F.6.2.b
    n1=0.955+0.006*psea; %Equation F.6.2a
    gammaw=(0.046+0.0019*psea+(3.98*n1)/((freq-22.235)^2+9.42*n1^2)*(1+((freq-22)/(freq+22))^2))*freq^2*psea*10^(-4); %Equation F.6.2

    if hsur<=2600
        psurr=psur+0.4+0.0003*hsur; %Equation F.5.1
    else 
        psurr=psur+5*exp(-hsur/1800); %Equation F.5.1
    end

    psur2=psurr; %Re-evaluate psur
    psear=psur2*exp(hsur/2000); %Equation F.6.2.b
    nr=0.955+0.006*psear; %Equation F.6.2a
    gammawr=(0.046+0.0019*psear+(3.98*nr)/((freq-22.235)^2+9.42*nr^2)*(1+((freq-22)/(freq+22))^2))*freq^2*psear*10^-4; %Equation F.6.2
end


%P.2001-2-2015: Attachment G4: Basic Transmission Loss
function [Lbe,foEs1,foEs2,Gamma1,Gamma2,LbEs1,LbEs2,Lp1r,Lp1t,Lp2r,Lp2t] = g_sporadic(p_time,phimn,phime,dn,freq,ae,phi1qe,phi1qn,phi3qe,phi3qn,thetat,thetar,dlt,dlr)
    %This method should not be considered reliable at low or high
    %geomagnetic latitudes, and it need not be calculated for a LoS path.
    
    %It should be noted that incidents of high signal strength due to this
    %phenomenon exhibit a very strong seasonal dependence.

    [LbEs1,foEs1,Gamma1,Lp1r,Lp1t] = g2_LbEs1(p_time,phimn,phime,dn,freq,ae,thetat,thetar,dlt,dlr);
    [LbEs2,foEs2,Gamma2,Lp2r,Lp2t] = g3_LbEs2(p_time,phi1qe,phi1qn,phi3qe,phi3qn,dn,freq,ae,thetat,thetar,dlt,dlr);
    
    Lbe=-10*log10(10^(-0.1*LbEs1)+10^(-0.1*LbEs2)); %Equation G.4.1c
    
    if LbEs1<(LbEs2-20)
        clear Lbe;
        Lbe=LbEs1; %Equation G.4.1a
    end
    
    if LbEs2<(LbEs1-20)
        clear Lbe;
        Lbe=LbEs2; %Equation G.4.1b
    end
end

%P.2001-2-2015: Attachment G1: Derivation of foEs
function [foEs] = g1_foEs(p_time,lat,lon)
    %This method should not be considered reliable at low or high
    %geomagnetic latitudes, and it need not be calculated for a LoS path.
    
    %It should be noted that incidents of high signal strength due to this
    %phenomenon exhibit a very strong seasonal dependence.

    FoEs50=csvread('nick_import_FoEs50.csv');
    FoEs10=csvread('nick_import_FoEs10.csv');
    FoEs01=csvread('nick_import_FoEs01.csv');
    FoEs0dot1=csvread('nick_import_FoEs0.1.csv');
    
    %Table G.1
    if p_time<1
       p1=0.1;
       p2=1;
       foEs1=get_txt_value(lat,lon,FoEs0dot1);
       foEs2=get_txt_value(lat,lon,FoEs01);
    end
    if p_time>=1 && p_time<=10
        p1=1;
        p2=10;
        foEs1=get_txt_value(lat,lon,FoEs01);
        foEs2=get_txt_value(lat,lon,FoEs10);
    end
    if p_time>10
        p1=10;
        p2=50;
        foEs1=get_txt_value(lat,lon,FoEs10);
        foEs2=get_txt_value(lat,lon,FoEs50);
    end
    
    foEs=foEs1+(foEs2-foEs1)*((log10(p_time/p1))/log10(p2/p1)); %Equation G.1.1 (MHz)
end


%P.2001-2-2015: Attachment G2: 1 Hop Propagation
function [LbEs1,foEs1,Gamma1,Lp1r,Lp1t] = g2_LbEs1(p_time,phimn,phime,dn,freq,ae,thetat,thetar,dlt,dlr)

    [foEs] = g1_foEs(p_time,phimn,phime); %Obtain foEs for the midpoint of the path.
    foEs1=foEs;
    Gamma1=((40/(1+(dn/130)+(dn/250)^2))+0.2*(dn/2600)^2)*((1000*freq)/foEs)^2+exp((dn-1660)/280); %Equation G.2.1

    hes=120; %hes is the height of the sporadic-E layer in km, set to 120km.

    l1=2*(ae^2+(ae+hes)^2-2*ae*(ae+hes)*cos(dn/(2*ae)))^(0.5); %Equation G.2.2

    Lbfs1=92.44+20*log10(freq)+20*log10(l1); %%Equation 3.11.1 %Equation G.2.3, free-space loss calculated for the slope distance.

    alpha1=dn/(2*ae); %Equation G.2.4a
    er1=0.5*pi-atan((ae*sin(alpha1))/(hes+ae*(1-cos(alpha1))))-alpha1; %Equation G.2.4

    d1t=0.001*thetat-er1; %Equation G.2.5
    d1r=0.001*thetar-er1; %Equation G.2.5
    
    if d1t>=0
        v1t=3.651*sqrt(1000*freq*dlt*((1-cos(d1t))/(cos(0.001*thetat)))); %Equation G.2.6a
    else
        v1t=-3.651*sqrt(1000*freq*dlt*((1-cos(d1t))/(cos(0.001*thetat)))); %Equation G.2.6b 
    end
    
    if d1r>=0
        v1r=3.651*sqrt(1000*freq*dlr*((1-cos(d1r))/(cos(0.001*thetar)))); %Equation G.2.6a
    else
        v1r=-3.651*sqrt(1000*freq*dlr*((1-cos(d1r))/(cos(0.001*thetar)))); %Equation G.2.6b 
    end
    
    Lp1t=knife_edge_Jv(v1t); %Equation G.2.7a
    Lp1r=knife_edge_Jv(v1r); %Equation G.2.7b
    
    LbEs1=Lbfs1+Gamma1+Lp1t+Lp1r; %Equation G.2.8
end

%P.2001-2-2015: Attachment G3: 2 Hop Propagation
function [LbEs2,foEs2h,Gamma2,Lp2r,Lp2t] = g3_LbEs2(p_time,phi1qe,phi1qn,phi3qe,phi3qn,dn,freq,ae,thetat,thetar,dlt,dlr)

    [foEs21] = g1_foEs(p_time,phi1qn,phi1qe); %Obtain foEs for the quarter point of the path
    [foEs22] = g1_foEs(p_time,phi3qn,phi3qe); %Obtain foEs for the quarter point of the path

    foEs2h=min([foEs21,foEs22]); 
    
    Gamma2=((40/(1+(dn/260)+(dn/500)^2))+0.2*(dn/5200)^2)*((1000*freq)/foEs2h)^2+exp((dn-3220)/560); %Equation G.3.1

    hes=120; %hes is the height of the sporadic-E layer in km, set to 120km.

    l2=4*(ae^2+(ae+hes)^2-2*ae*(ae+hes)*cos(dn/(4*ae)))^(0.5); %Equation G.3.2

    Lbfs2=92.44+20*log10(freq)+20*log10(l2); %%Equation 3.11.1 %Equation G.3.3, free-space loss calculated for the slope distance.

    alpha2=dn/(4*ae); %Equation G.3.4a
    er2=0.5*pi-atan((ae*sin(alpha2))/(hes+ae*(1-cos(alpha2))))-alpha2; %Equation G.3.4

    d2t=0.001*thetat-er2; %Equation G.3.5
    d2r=0.001*thetar-er2; %Equation G.3.5
    
    if d2t>=0
        v2t=3.651*sqrt(1000*freq*dlt*((1-cos(d2t))/(cos(0.001*thetat)))); %Equation G.3.6a
    else
        v2t=-3.651*sqrt(1000*freq*dlt*((1-cos(d2t))/(cos(0.001*thetat)))); %Equation G.3.6b 
    end
    
    if d2r>=0
        v2r=3.651*sqrt(1000*freq*dlr*((1-cos(d2r))/(cos(0.001*thetar)))); %Equation G.3.6a
    else
        v2r=-3.651*sqrt(1000*freq*dlr*((1-cos(d2r))/(cos(0.001*thetar)))); %Equation G.3.6b 
    end
    
    Lp2t=knife_edge_Jv(v2t); %Equation G.3.7a
    Lp2r=knife_edge_Jv(v2r); %Equation G.3.7b
    
    LbEs2=Lbfs2+Gamma2+Lp2t+Lp2r; %Equation G.3.8
end

%P.2001-2-2015: Attachement I: Iterative Procedure to invert a cumulative distribution function.
function [Aiter] = I_Aiter_submodel1(Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1,q_time,Q0ca)

    %Section I.2
    %Stage 1, setting the search range
    Ainit=10; %10dB
    Ahigh=Ainit/2; %Equation I.2.1
    Alow=-Ainit/2; %Equation I.2.2
    Astep=Ainit; %Equation I.2.3

    qhigh=Qiter_sm1(Ahigh,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1); %Equation I.2.4a
    qlow=Qiter_sm1(Alow,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1); %Equation I.2.4b

    q=q_time;

    %Stage 1 
    for i=1:1:10 %No more than 10 loops
        if q<qhigh
            Alow=Ahigh;
            qlow=qhigh;
            Astep=2*Astep;
            Ahigh=Ahigh+Astep;
            qhigh=Qiter_sm1(Ahigh,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1);
        end

        if q>qlow
            Ahigh=Alow;
            qhigh=qlow;
            Astep=2*Astep;
            Alow=Alow-Astep;
            qlow=Qiter_sm1(Alow,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1);
        end

        if q>=qhigh && q<=qlow
           break; %Proceed to Stage 2
        end
        
    end

    %Stage 2: Binary Search
    
    Atry=0.5*(Alow+Ahigh); %Equation I.2.5
    qtry=Qiter_sm1(Atry,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1);  %Equation I.2.6

    Aacc=0.01;
    niter=ceil(3.32*log10(Astep/Aacc));

    for i=1:1:niter
        if qtry<q
            Ahigh=Atry;
        else
            Alow=Atry;
        end
        Atry=0.5*(Alow+Ahigh); %Equation I.2.5
        qtry=Qiter_sm1(Atry,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1);  %Equation I.2.6
    end

    Aiter=Atry; 
end

%P.2001-2-2015: Qiter(q) for Submodel 1
function [varout] = Qiter_sm1(A,Q0ca,Q0ra,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1) 

    [Qcaf] = b4_Qcaf(Q0ca,A); %B.4: Calculate Qcaf(A)

    [Qrain]=c3_Qrain(A,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1); %C.3: Calculate Qrain(A)

    varout=Qrain*(Q0ra/100)+Qcaf*(1-(Q0ra/100)); %Qiter(A):Equation 4.1.3
    
end

%P.2001-2-2015: Attachement I: Iterative Procedure to invert a cumulative distribution function.
%It seems like there should be a way to do this with Matlab
function [Aiter] = I_Aiter_submodel3(Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1,q_time)

    %Section I.2
    %Stage 1, setting the search range
    Ainit=10; %10dB
    Ahigh=Ainit/2; %Equation I.2.1
    Alow=-Ainit/2; %Equation I.2.2
    Astep=Ainit; %Equation I.2.3

    qhigh=Qiter_sm3(Ahigh,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1); %Equation I.2.4a
    qlow=Qiter_sm3(Alow,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1); %Equation I.2.4b
    
    q=q_time;

    %Stage 1 
    for i=1:1:10 %No more than 10 loops
        if q<qhigh
            Alow=Ahigh;
            qlow=qhigh;
            Astep=2*Astep;
            Ahigh=Ahigh+Astep;
            qhigh=Qiter_sm3(Ahigh,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1);
        end

        if q>qlow
            Ahigh=Alow;
            qhigh=qlow;
            Astep=2*Astep;
            Alow=Alow-Astep;
            qlow=Qiter_sm3(Alow,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1);
        end

        if q>=qhigh && q<=qlow
           break; %Proceed to Stage 2
        end 
    end

    %Stage 2: Binary Search
    
    Atry=0.5*(Alow+Ahigh); %Equation I.2.5
    qtry=Qiter_sm3(Atry,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1);  %Equation I.2.6

    Aacc=0.01;
    niter=ceil(3.32*log10(Astep/Aacc));

    for i=1:1:niter
        if qtry<q
            Ahigh=Atry;
        else
            Alow=Atry;
        end
        Atry=0.5*(Alow+Ahigh); %Equation I.2.5
        qtry=Qiter_sm3(Atry,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1);  %Equation I.2.6
    end

    Aiter=Atry; 
end

%P.2001-2-2015: Qiter(q) for Submodel 3
function [varout] = Qiter_sm3(A,Q0ra,Pr6,Pm,Gm,hRtop,hrainlo,dr,kmod,alphamod,a1,b1,c1) 
    Qcaftropo=b5_Qcaftropo(A);

    Qrain=c3_Qrain(A,hrainlo,hRtop,Pm,Gm,Pr6,dr,kmod,alphamod,a1,b1,c1); %C.3: Calculate Qrain(A)

    varout=Qrain*(Q0ra/100)+Qcaftropo*(1-(Q0ra/100)); %Qiter(A):Equation 4.1.3    
end

%P.2001-2-2015: 3.12: Knife-Edge Diffraction Loss
function [Jv] = knife_edge_Jv(v)
    if v>(-0.78)
    	Jv=6.9+(20*log10(sqrt(((v-0.1)^2)+1)+v-0.1)); %Equation 3.12.1a 
    else
    	Jv=0; %Equation 3.12.1b
    end
end

%P.2001-2-2015: Find Value in txt data
function [value] = get_txt_value(lat,lon,temp) 
    %Based on the size of the data, find the 4 nearest points and Bilinear Interp
    %This is only valid with nick_import_"data".csv
    data=flipud(temp);
    [y3, x3]=size(data);

    if x3==720 %For Tropo Data, take nearest point, no bi-linear interp
       temp_lon=lon+180; %Make Longitude from 0-360  
       temp_lat=lat+90;  %Make latitude from 0-180
       spacing_lat=180/(y3);
       spacing_lon=360/(x3);

       %find two nearest lats and lons
        series_lon=1:1:x3;
        series_lat=1:1:y3;

        %Makes it larger
        small_lat=temp_lat*(y3)/180;
        small_lon=temp_lon*(x3)/360;

        %Find 2 nearest lats/longs to small lat/lon
        [lon_val, lon_idx]=sort(abs(small_lon+0.25-series_lon));
        [lat_val, lat_idx]=sort(abs(small_lat+0.25-series_lat));

        lon2_idx=sort(lon_idx(1));
        lat2_idx=sort(lat_idx(1));

        value=data(series_lat(lat2_idx),series_lon(lon2_idx));

    else %All other data files
        spacing_lat=180/(y3-1);
        spacing_lon=360/(x3-1);


        temp_lat=lat+90;  %Make latitude from 0-180
        if lon<0 %Make longitude from 0-360
            temp_lon=lon+360; 
        else
            temp_lon=lon;
        end

        %find two nearest lats and lons
        series_lon=1:1:x3;
        series_lat=1:1:y3;

        %Makes it smaller
        small_lat=temp_lat*(y3-1)/180+1;
        small_lon=temp_lon*(x3-1)/360+1;

        %Find 2 nearest lats/longs to small lat/lon
        [lon_val, lon_idx]=sort(abs(small_lon-series_lon));
        [lat_val, lat_idx]=sort(abs(small_lat-series_lat));

        lon2_idx=sort(lon_idx(1:2));
        lat2_idx=sort(lat_idx(1:2));

        Q11=data(series_lat(lat2_idx(1)),series_lon(lon2_idx(1)));
        Q12=data(series_lat(lat2_idx(2)),series_lon(lon2_idx(1)));
        Q21=data(series_lat(lat2_idx(1)),series_lon(lon2_idx(2)));
        Q22=data(series_lat(lat2_idx(2)),series_lon(lon2_idx(2)));

        x=small_lon;
        y=small_lat;

        x1=series_lon(lon2_idx(1));
        x2=series_lon(lon2_idx(2));
        y1=series_lat(lat2_idx(1));
        y2=series_lat(lat2_idx(2));

        value=((((x2-x)*(y2-y))/((x2-x1)*(y2-y1)))*Q11)+((((x-x1)*(y2-y))/((x2-x1)*(y2-y1)))*Q21)+((((x2-x)*(y-y1))/((x2-x1)*(y2-y1)))*Q12)+((((x-x1)*(y-y1))/((x2-x1)*(y2-y1)))*Q22); %This is the return value
    end
end


%P.2001: D.1: Find the Distance for dlm/dtm
function [dtm,dlm,dct,dcr,omega] = find_ITU_dist(rx_pos,tx_pos,h,dn)

    TropoClim=csvread('nick_import_TropoClim.csv'); 

    ellipsoid=[];
    mid_npts_rx=length(h);% This is based upon the number of points in the terrain model (h).
    [r_lat,r_lon]=track2(rx_pos(1),rx_pos(2),tx_pos(1),tx_pos(2),ellipsoid,'degrees',mid_npts_rx);
    for i=1:1:length(r_lat)
        tropo_value_rx(i)=get_txt_value(r_lat(i),r_lon(i),TropoClim); 
        %Gets the Climate of each point starting at the rx
    end
        
    clear consec_dist_rx;
    marker1=1;
    ind1=1;
    ind2=0;
    for i=1:1:length(tropo_value_rx)
        if tropo_value_rx(i)~=0
            ind2=i;
        else
            if ind2==0
                ind2=i;
            else
                consec_dist_rx(marker1)=deg2km(distance(r_lat(ind1),r_lon(ind1),r_lat(ind2),r_lon(ind2)));
                marker1=marker1+1;
                ind1=i;
                if i<length(tropo_value_rx)
                    ind1=i+1;
                end
            end
        end
    end
    consec_dist_rx(marker1)=deg2km(distance(r_lat(ind1),r_lon(ind1),r_lat(ind2),r_lon(ind2)));
    
    consec_len_rx=max(consec_dist_rx);
    
    dcr=consec_dist_rx(1); %The r_lat/lon starts at the receiver
    
    
    mid_npts_tx=length(h);% This is based upon the number of points in the terrain model (h).
    [t_lat,t_lon]=track2(tx_pos(1),tx_pos(2),rx_pos(1),rx_pos(2),ellipsoid,'degrees',mid_npts_rx);
    for i=1:1:length(t_lat)
        tropo_value_tx(i)=get_txt_value(t_lat(i),t_lon(i),TropoClim); 
        %Gets the Climate of each point starting at the rx
    end
    
    clear consec_dist_tx;
    marker1=1;
    ind1=1;
    ind2=0;
    for i=1:1:length(tropo_value_tx)
        if tropo_value_tx(i)~=0
            ind2=i;
        else
            if ind2==0
                ind2=i;
            else
                consec_dist_tx(marker1)=deg2km(distance(t_lat(ind1),t_lon(ind1),t_lat(ind2),t_lon(ind2)));
                marker1=marker1+1;
                ind1=i;
                if i<length(tropo_value_tx)
                    ind1=i+1;
                end
            end
        end
    end
    consec_dist_tx(marker1)=deg2km(distance(t_lat(ind1),t_lon(ind1),t_lat(ind2),t_lon(ind2)));
    
    consec_len_tx=max(consec_dist_tx);
    
     dct=consec_dist_tx(1); %The t_lat/lon starts at the receiver

     %Need to refine how to define dtm, dlm based on h.
     
     %Find the first inland point near the sea, and go inland, checking the
     %elevation to see when it become greater than 100m, and then see how long that path if,
     %checking that it does not exceed 50km.
     
     
     
    dtm=consec_len_tx;  %Coast and inland distance.
    dlm=consec_len_tx;  %Inland A2, greater than 100m, not to exceed 50km inland.
    
    omega=round(1-(consec_len_tx/dn),2);
end



















