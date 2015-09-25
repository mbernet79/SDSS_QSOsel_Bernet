function[result_ellipse,distance_caxis,flag_special,min_index]=stellar_locus_ugri_bu(umag,gmag,rmag,imag,erru,errg,errr,erri,N_sigma,imag_cut)
%Function to determine whether an object with given u,g,r,i magnitudes
%is located outside of the stellar locus. 
%Implementation is according to algorithm outlined in Richards G. T. et al., 2002, AJ, 123, 2945
%Parametrization of stellar locus is explained in Newberg, H.J. & Yanny, B.
%1997, ApJS, 113,89. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flag_special=0: everything is o.k.
%flag_special=1: outlier of the half ellipsoids
%flag_special=2: outside of magnitude cut

flag_special=0;
result_ellipse=0;

%Check if outside imag cut

if(imag >= imag_cut)      
    result_ellipse=0;
    distance_caxis=0.0;
    flag_special=2;
    min_index=0;
    return
else

%Parameters for the endpoints
    
k_blue=-0.05;
a_kblue=0.2; 
if(N_sigma < 4.0)
    a_kblue=0.1;
end
k_red=100.0;
a_kred=0.0;

%position vector of QSO in color-color-color space

r_QSO(1,1)=umag-gmag;         
r_QSO(2,1)=gmag-rmag;
r_QSO(3,1)=rmag-imag;


%Calculate matrix V_ab as defined in Richards et al. (2002):
%V_ab=Q^T*S*Q, where S is the covariance matrix of (u-g,g-r,r-i), in the 3D
%color space and V_ab is the 2D covariance matrix in the (l,m) plane (major axis, minor axis of cylinders). 
%Q is the corresponding transformation matrix.

S(1:3,1:3)=0.0;
S(1,1)=erru^2+errg^2;
S(2,2)=errg^2+errr^2;
S(3,3)=errr^2+erri^2;

S(1,2)=-errg^2;
S(2,3)=-errr^2;
S(2,1)=S(1,2);
S(3,2)=S(2,3);

%vector containing the positions of the locus points (17 points)
position_lp(:,1)=[0.855;1.002;1.136;1.262;1.382;1.499;1.609;1.712;1.808;1.904;2.007;2.117;2.234;2.361;2.478;2.518;2.510]; %u-g
position_lp(:,2)=[0.259;0.344;0.410;0.466;0.517;0.565;0.611;0.655;0.700;0.748;0.802;0.866;0.945;1.047;1.191;1.327;1.355]; %g-r
position_lp(:,3)=[0.094;0.126;0.150;0.170;0.189;0.207;0.223;0.238;0.255;0.273;0.293;0.317;0.351;0.403;0.502;0.707;1.068]; %r-i

%k-vectors of the locus (17 points)

kv(:,1)=[0.851;0.868;0.893;0.907;0.915;0.915;0.912;0.904;0.888;0.878;0.860;0.827;0.773;0.646;0.355;0.053;-0.022]; %kx  before k(17,1) = -0.022
kv(:,2)=[0.492;0.467;0.422;0.396;0.379;0.379;0.387;0.402;0.430;0.449;0.478;0.521;0.573;0.650;0.634;0.278;0.076];  %ky  before k(17,2) = 0.076
kv(:,3)=[0.182;0.172;0.154;0.145;0.140;0.135;0.132;0.146;0.165;0.166;0.178;0.213;0.271;0.400;0.688;0.959;0.997];  %kz  before k(17,3) = 0.997

%vector cotaining locus ellipse parameters; (17 points)

le(:,1)=[0.282;0.247;0.221;0.219;0.216;0.217;0.224;0.227;0.233;0.248;0.266;0.278;0.309;0.382;0.463;0.484;0.569]; %major axis, l
le(:,2)=[0.135;0.129;0.124;0.126;0.125;0.129;0.131;0.127;0.132;0.129;0.134;0.136;0.136;0.145;0.156;0.180;0.212]; %minor axis, m
le(:,3)=[-1.165;-1.147;-1.075;-1.026;-0.977;-0.983;-0.986;-0.989;-1.040;-1.002;-1.017;-1.023;-1.033;-1.051;-1.108;-1.244;-1.669]; %position angle ellipse %before -1.669

%Finding locus point with minimum distance

for i=1:length(position_lp(:,1))
    
    distance(i,1)=sqrt((r_QSO(1,1)-position_lp(i,1))^2 + (r_QSO(2,1)-position_lp(i,2))^2 + (r_QSO(3,1)-position_lp(i,3))^2);
end

[min_dist,min_index]=min(distance(:,1));

sqrt_kx_ky=sqrt(kv(min_index,1)^2+kv(min_index,2)^2);

%Calculating the entry of Q

Q(1,1)= (sin(le(min_index,3))*kv(min_index,2)- cos(le(min_index,3))*kv(min_index,1)*kv(min_index,3))/sqrt_kx_ky;
Q(2,1)= (-sin(le(min_index,3))*kv(min_index,1) - cos(le(min_index,3))*kv(min_index,2)*kv(min_index,3))/sqrt_kx_ky;
Q(3,1)= cos(le(min_index,3))*sqrt_kx_ky;
Q(1,2)= (cos(le(min_index,3))*kv(min_index,2) + sin(le(min_index,3))*kv(min_index,1)*kv(min_index,3))/sqrt_kx_ky;
Q(2,2)= (-cos(le(min_index,3))*kv(min_index,1) + sin(le(min_index,3))*kv(min_index,2)*kv(min_index,3))/sqrt_kx_ky;
Q(3,2)= -sin(le(min_index,3))*sqrt_kx_ky;

Q_T=transpose(Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_ab=Q_T*S*Q;

amajor_pow2=N_sigma^2*(V_ab(1,1)+V_ab(2,2)+sqrt((V_ab(1,1)-V_ab(2,2))^2+4*V_ab(1,2)^2))/2.0;
amajor=sqrt(amajor_pow2);

aminor_pow2=N_sigma^2*(V_ab(1,1)+V_ab(2,2)-sqrt((V_ab(1,1)-V_ab(2,2))^2+4*V_ab(1,2)^2))/2.0;
aminor=sqrt(aminor_pow2);
if(V_ab(1,2)==0 & V_ab(1,1)==0 & V_ab(2,2)==0)
    theta_err=0.0;
end
if(V_ab(1,2)==0 & V_ab(1,1)>V_ab(2,2))
    theta_err=0.0;
end
if(V_ab(1,2)==0 & V_ab(2,2)>V_ab(1,1))
    theta_err=pi/2;
end
if(V_ab(1,2)~=0)
    tan_theta=(-(V_ab(1,1)-V_ab(2,2))+sqrt((V_ab(1,1)-V_ab(2,2))^2+4*V_ab(1,2)^2))/(2*V_ab(1,2));
    theta_err=atan(tan_theta);
end
if(V_ab(1,2)==0 & V_ab(1,1)==0 & V_ab(2,2)==0.0)
    theta_err=0.0;
end
    
%Here are the parameters of the convolved ellipse calculated

alpha=aminor^2*cos(theta_err)^2+amajor^2*sin(theta_err)^2+le(min_index,2)^2;

beta=sin(theta_err)*cos(theta_err)*(amajor^2-aminor^2);

gamma=aminor^2*sin(theta_err)^2+amajor^2*cos(theta_err)^2+le(min_index,1)^2;

a_r=(alpha+gamma)/2;

b_r=sqrt((alpha+gamma)^2-4*(alpha*gamma-beta^2))/2;

a_total_major=sqrt(a_r+b_r);

a_total_minor=sqrt(a_r-b_r);

% if(beta==0 & alpha==gamma) % I suspect these are the wrong conditions and
%     theta_err_total=0;     % implemented instead the condition below. 
% end
% if(beta==0 & alpha~=gamma)
%     theta_err_total=pi/2;
% end
if(beta==0)
   theta_err_total=theta_err;   
end
if(beta~=0 & theta_err >=0)
    theta_err_total=atan((alpha-gamma+sqrt((alpha-gamma)^2+4*beta^2))/(2*beta));
end
if(beta~=0 & theta_err < 0)
    theta_err_total=-atan((alpha-gamma+sqrt((alpha-gamma)^2+4*beta^2))/(2*beta));
end

theta_err_total_final=le(min_index,3)+theta_err_total;    % theta_err_total_final is defined with respect to i vector

if(theta_err_total_final < -pi/2)                         % Be sure that theta_err_total_final is within -pi/2 and +pi/2 
    
    theta_err_total_final=theta_err_total_final +pi;
end
if(theta_err_total_final > pi/2)
    
    theta_err_total_final=theta_err_total_final -pi;
end

%Calculate increase of end ellipsoids at blue and red end of locus

Q_k_blue(1,1)=kv(1,1);
Q_k_blue(2,1)=kv(1,2);
Q_k_blue(3,1)=kv(1,3);

Var_k_blue=(Q_k_blue.')*S*Q_k_blue;

Q_k_red(1,1)=kv(17,1);
Q_k_red(2,1)=kv(17,2);
Q_k_red(3,1)=kv(17,3);

Var_k_red=(Q_k_red.')*S*Q_k_red;

%Calculating distance from axis of cylinder and comparing with extension of
%cylinder

k_locus_h=kv(min_index,:);    %Redefining some things for easier calculation
k_locus=k_locus_h.';

r_P(1,1)=position_lp(min_index,1);  % locus point which is closest to r_QSO
r_P(2,1)=position_lp(min_index,2);
r_P(3,1)=position_lp(min_index,3);

%third vector according to Newberg (not needed)
temp=sqrt(k_locus(1,1).*k_locus(1,1)+k_locus(2,1).*k_locus(2,1));  % writing out the i_vector which is a unit vector in the plane perpendicular to k.
i_vector(1,1)=-k_locus(1,1).*k_locus(3,1)./temp;
i_vector(2,1)=-k_locus(2,1).*k_locus(3,1)./temp;
i_vector(3,1)=temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=(k_locus.')*(r_QSO-r_P);
r_A=r_P+lambda.*k_locus;            % point on cylinder axis which is closest to r_QSO
distance_caxis= sqrt(((r_QSO-r_A).')*(r_QSO-r_A));   % distance of r_QSO to cylinder axis
Theta_QSO=acos(((r_QSO-r_A).')*i_vector/distance_caxis);

%Calculate radius of ellipse at angle theta_err_total_final-Theta_QSO

eccentricity=sqrt(a_total_major.*a_total_major-a_total_minor.*a_total_minor)./a_total_major;
r_ellipse=a_total_minor./sqrt(1-eccentricity.^2.*cos(theta_err_total_final-Theta_QSO).^2);

%Calculate if |r_A-r_P| for min_index=1 or 17 is larger than a_kblue or a_kred
%and calculate if r_QSO is outside the half ellipsoid

distance_r_A=0.0;
if((min_index==1 & lambda < k_blue)|(min_index==17 & lambda > k_red))
    
    flag_special=1;
    if(min_index==1)
        a_ktot=sqrt(a_kblue.*a_kblue+N_sigma.*N_sigma.*Var_k_blue);
        k_tot=k_blue;
        r_End=r_P+k_blue.*k_locus;
    end
    if(min_index==17)
        a_ktot=sqrt(a_kred.*a_kred+N_sigma.*N_sigma.*Var_k_red);
        k_tot=k_red;
        r_End=r_P+k_red.*k_locus;
    end
    distance_r_A= sqrt(((r_End-r_A).')*(r_End-r_A)); 
    angle_u=atan2(distance_caxis,distance_r_A);
    angle_v=theta_err_total_final-Theta_QSO;
   
    if(a_ktot~=0)
        
    a_total_major_end=a_total_major*sqrt(1-distance_r_A.^2/a_ktot.^2);  %These are the major and minor semi axes of the ellipses at distance z from the ellipsoid
    a_total_minor_end=a_total_minor*sqrt(1-distance_r_A.^2/a_ktot.^2);
    else
        a_total_major_end=0;
        a_total_minor_end=0;
    end
    eccentricity=sqrt(a_total_major_end.*a_total_major_end-a_total_minor_end.*a_total_minor_end)./a_total_major_end;
    r_ellipse=a_total_minor_end./sqrt(1-eccentricity.^2.*cos(theta_err_total_final-Theta_QSO).^2);
       
    if(distance_caxis > r_ellipse)
           
           result_ellipse=1;
    end
end

%Check if r_QSO is outside the ellipse

if(distance_caxis > r_ellipse & flag_special==0) 
    result_ellipse=1;
    
end



end




