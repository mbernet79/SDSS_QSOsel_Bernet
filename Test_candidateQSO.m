%This script tests whether an object with given u,g,r,i,z magnitudes and associated
%errors is selected by the SDSS QSO selection algorithm as outlined in
%Richards G. T. et al., 2002, AJ, 123, 2945
%
%Please check the readme.txt file for the applicability of the algorithm and its
%restrictions
%
%Please cite Murphy, M.T., Bernet, M.L., 2015, MNRAS, submitted, if you use this code
%Comments or reports of errors are welcome. Email: bernet.martin@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umag=18.0;  % e.g. 18.0
gmag=18.0;   
rmag=18.0;  
imag=18.0;  
zmag=18.0;  

err_umag=0.02;      %e.g. 0.02
err_gmag=0.02;  
err_imag=0.02;  
err_rmag=0.02;  
err_zmag=0.02;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[result(:,1),temp1,temp2,result(:,2)]=stellar_locus_ugri(umag,gmag,rmag,imag,err_umag,err_gmag,err_imag,err_zmag,4.0,19.1);
     
[result(:,3),temp3,temp4,result(:,4)]=stellar_locus_griz(umag,gmag,rmag,imag,zmag,err_gmag,err_rmag,err_imag,err_zmag,4.0);
%flag low_z    high_z
[result(:,5),result(:,6),result(:,7)]=flagging_incl_region(umag,gmag,rmag,imag,zmag,err_umag,err_gmag,err_rmag,err_imag);
     
[result(:,8),result(:,9)]=flagging_excl_region(umag,gmag,rmag,imag,zmag,err_umag,err_gmag,err_rmag,err_imag,err_zmag);
     
