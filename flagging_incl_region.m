function[flag_incl_lowz,flag_incl_highz,code_incl]=flagging_incl_region(umag,gmag,rmag,imag,zmag,erru,errg,errr,erri)
%Function to check whether an object with given u,g,r,i,z magnitudes is in
%an inclusion region. If yes, either flag_incl_lowz or flag_incl_highz is
%set to 1.


flag_incl_lowz=0;
flag_incl_highz=0;
code_incl=0;


% mid z inclusion region
% Check whether you want to include only every tenth object (uncomment below)
if(umag-gmag > 0.6 && umag-gmag < 1.5)                
    if(gmag-rmag > 0.0 && gmag-rmag < 0.2)
        if(rmag-imag > -0.1 && rmag -imag < 0.4)
            if(imag-zmag > -0.1 && imag -zmag < 0.4)
                [flag_midz,temp_x1,temp_x2,temp_x3]=stellar_locus_ugri(umag,gmag,rmag,imag,erru,errg,errr,erri,2.0,19.1);
                if(flag_midz==1)% && randi(10,1)==7)    %only selects every tenth object which fullfills this criteria
                    flag_incl_lowz=1;
                    code_incl=1;
                end
            end
        end
    end
end

%UVX inclusion region
if(umag > 0.0 && gmag > 0.0 && errg < 0.1 && errr < 0.1 && umag - gmag < 0.6) 
    flag_incl_lowz=1;
    code_incl=2;
end

% z > 3.6 inclusion region
if(erri < 0.2)                 
    if(umag-gmag > 1.5 | umag > 20.6)
        if(gmag -rmag > 0.7)
            if(gmag -rmag > 2.1 | rmag-imag < 0.44*(gmag-rmag)-0.358)
                if(imag-zmag > -1.0 && imag -zmag < 0.25)
                    flag_incl_highz=1;
                    code_incl=3;
                end
            end
        end
    end
end

% z > 4.5 inclusion region
if(erri < 0.2)                
    if(umag > 21.5 && gmag > 21.0)
        if(rmag-imag > 0.6)
            if(imag-zmag > -1.0 && imag -zmag < 0.52*(rmag-imag)-0.412)
                flag_incl_highz=1;
                code_incl=4;
            end
        end
    end
end

%ugr red outlier inclusion region

[flag_ugr,temp_x1,temp_x2,temp_x3]=stellar_locus_ugri(umag,gmag,rmag,imag,erru,errg,errr,erri,4.0,20.2);
                
if((umag > 20.6 && umag-gmag > 1.5 && gmag-rmag < 1.2 && rmag-imag < 0.3 && imag-zmag > -1.0 && (gmag-rmag < 0.44*(umag-gmag)-0.56)) | (flag_ugr==1 && umag > 0.0 && gmag > 0.0 && errg < 0.2 && errr < 0.2 && umag - gmag > 1.5)) 
    flag_incl_highz=1;
    code_incl=5;
end


end

