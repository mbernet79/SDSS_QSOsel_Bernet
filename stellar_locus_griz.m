function[result_ellipse,distance_caxis,flag_special,min_index]=stellar_locus_griz(umag,gmag,rmag,imag,zmag,errg,errr,erri,errz,N_sigma)
%Function to determine whether an object with given g,r,i,z magnitudes
%is located outside of the stellar locus. 
%Implementation is according to algorithm outlined in Richards G. T. et al., 2002, AJ, 123, 2945
%Parametrization of stellar locus is explained in Newberg, H.J. & Yanny, B.
%1997, ApJS, 113,89. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%flag_special=0: everything is o.k.
%flag_special=1: outlier of the half ellipsoids
%flag_special=2: outside of magnitude cut
%flag_special=3: griz outlier, but not selected because potential lowz qso

flag_special=0;
result_ellipse=0;

%Parameters for the endpoints
k_blue=-0.3;
a_kblue=0.5; 
k_red=100.0;
a_kred=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check if outside imag cut
if(imag >= 20.2)      
    result_ellipse=0;
    distance_caxis=0.0;
    flag_special=2;
    min_index=0;
    return
else

    
%position vector of QSO in color-color-color space

r_QSO(1,1)=gmag-rmag;         
r_QSO(2,1)=rmag-imag;
r_QSO(3,1)=imag-zmag;

%Calculate matrix V_ab as defined in Richards et al. (2002):
%V_ab=Q^T*S*Q, where S is the covariance matrix of (u-g,g-r,r-i), in the 3D
%color space and V_ab is the 2D covariance matrix in the (l,m) plane (major axis, minor axis of cylinders). 
%Q is the corresponding transformation matrix.

S(1:3,1:3)=0.0;
S(1,1)=errg^2+errr^2;
S(2,2)=errr^2+erri^2;
S(3,3)=erri^2+errz^2;

S(1,2)=-errr^2;
S(2,3)=-erri^2;
S(2,1)=S(1,2);
S(3,2)=S(2,3);

%vector containing the positions of the locus points (23 points)
position_lp(:,1)=[0.204;0.304;0.382;0.454;0.525;0.594;0.659;0.723;0.787;0.853;0.922;0.991;1.063;1.132;1.202;1.262;1.313;1.343;1.355;1.352;1.347;1.350;1.361]; %g-r
position_lp(:,2)=[0.071;0.110;0.137;0.166;0.194;0.219;0.242;0.265;0.288;0.313;0.341;0.371;0.409;0.454;0.507;0.569;0.651;0.754;0.874;0.996;1.116;1.240;1.385]; %r-i
position_lp(:,3)=[0.003;0.027;0.044;0.066;0.087;0.105;0.123;0.140;0.155;0.171;0.188;0.206;0.227;0.251;0.280;0.314;0.356;0.408;0.465;0.525;0.583;0.646;0.729]; %i-z

%k-vectors of the locus (23 points)

kv(:,1)=[0.911;0.916;0.910;0.895;0.905;0.913;0.911;0.915;0.916;0.906;0.897;0.876;0.832;0.778;0.704;0.566;0.362;0.168;0.035;-0.031;-0.008;0.047;0.067]; %kx
kv(:,2)=[0.351;0.339;0.340;0.356;0.342;0.325;0.330;0.332;0.335;0.360;0.380;0.420;0.485;0.551;0.623;0.729;0.832;0.885;0.900;0.899;0.895;0.879;0.868];  %ky
kv(:,3)=[0.218;0.213;0.237;0.268;0.253;0.246;0.246;0.231;0.220;0.222;0.227;0.237;0.267;0.301;0.342;0.386;0.420;0.434;0.435;0.438;0.446;0.475;0.493];  %kz

%vector cotaining locus ellipse parameters; (23 points)

le(:,1)=[0.207;0.165;0.154;0.159;0.164;0.162;0.151;0.150;0.153;0.157;0.160;0.163;0.171;0.175;0.178;0.185;0.193;0.213;0.246;0.250;0.265;0.246;0.300]; %major axis, l
le(:,2)=[0.146;0.126;0.128;0.134;0.133;0.133;0.133;0.127;0.128;0.124;0.125;0.123;0.123;0.125;0.127;0.135;0.129;0.131;0.137;0.135;0.133;0.121;0.139]; %minor axis, m
le(:,3)=[0.067;-2.907;-2.990;-0.029;-0.194;-0.315;-0.610;-0.858;-0.935;-0.917;-0.921;-0.898;-0.949;-1.033;-1.127;-1.323;-1.423;-1.554;-1.628;-1.667;-1.647;-1.652;-1.530]; %position angle ellipse

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%if(beta==0 & alpha==gamma) % I suspect these are the wrong conditions and
%    theta_err_total=0;     % implemented instead the condition below. 
%end
%if(beta==0 & alpha~=gamma) 
%    theta_err_total=pi/2;
%end

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

Q_k_red(1,1)=kv(23,1);
Q_k_red(2,1)=kv(23,2);
Q_k_red(3,1)=kv(23,3);

Var_k_red=(Q_k_red.')*S*Q_k_red;

%Calculating distance from axis of cylinder and comparing with extension of
%cylinder

k_locus_h=kv(min_index,:);    %Redefining some things for easier calculation
k_locus=k_locus_h.';

r_P(1,1)=position_lp(min_index,1);  % locus point which is closest to r_QSO
r_P(2,1)=position_lp(min_index,2);
r_P(3,1)=position_lp(min_index,3);

temp=sqrt(k_locus(1,1).*k_locus(1,1)+k_locus(2,1).*k_locus(2,1));  % writing out the i_vector which is a unit vector in the plane perpendicular to k.
%third vector according to Newberg (not needed)
i_vector(1,1)=-k_locus(1,1).*k_locus(3,1)./temp;
i_vector(2,1)=-k_locus(2,1).*k_locus(3,1)./temp;
i_vector(3,1)=temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=(k_locus.')*(r_QSO-r_P);

r_A=r_P+lambda.*k_locus;            % point on cylinder axis which is closest to r_QSO
distance_caxis= sqrt(((r_QSO-r_A).')*(r_QSO-r_A));   % distance of r_QSO to cylinder axis
Theta_QSO=acos(((r_QSO-r_A).')*i_vector/distance_caxis);  % angle with respect to r_QSO and i_vector

%Calculate radius of ellipse at angle theta_err_total_final-Theta_QSO

eccentricity=sqrt(a_total_major.*a_total_major-a_total_minor.*a_total_minor)./a_total_major;
r_ellipse=a_total_minor./sqrt(1-eccentricity.^2.*cos(theta_err_total_final-Theta_QSO).^2);


%Calculate if |r_A-r_P| for min_index=1 or 23 is larger than a_kblue or a_kred
%and calculate if r_QSO is outside the half ellipsoid

distance_r_A=0.0;
if((min_index==1 & lambda < k_blue)|(min_index==23 & lambda > k_red))
    
    flag_special=1;
    if(min_index==1)
        a_ktot=sqrt(a_kblue.*a_kblue+N_sigma.*N_sigma.*Var_k_blue);
        k_tot=k_blue;
        r_End=r_P+k_blue.*k_locus;
    end
    if(min_index==23)
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
        
        % This cut is applied to avoid that high z selection is dominated by low z QSOs   
           
        if(gmag-rmag < 1.0 & umag-gmag >= 0.8)  
    
            if(imag >=19.1 | umag-gmag < 2.5)
        
            flag_special=3;    
            result_ellipse=0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end


%Check if r_QSO is outside the ellipse

if(distance_caxis > r_ellipse & flag_special==0) 

    result_ellipse=1;    
    % This cut is applied to avoid that high z selection is dominated by
    % lowz QSOs
    if(gmag-rmag < 1.0 & umag-gmag >= 0.8)  
        if(imag >=19.1 | umag-gmag < 2.5)
        
            flag_special=3;    
            result_ellipse=0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end





