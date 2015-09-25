function[flag_excl,code_excl]=flagging_excl_region(umag,gmag,rmag,imag,zmag,erru,errg,errr,erri,errz)
%Function to check whether an object with given u,g,r,i magnitudes and associated errors
%lies in an exclusion region


flag_excl=0;
code_excl=0;    

%White dwarfs exclusion region

if(umag-gmag > -0.8 & umag-gmag < 0.7)    
    if(gmag-rmag > -0.8 & gmag-rmag < -0.1)
        if(rmag-imag > -0.6 & rmag-imag < -0.1)
            if(imag-zmag > -1.0 & imag-zmag < -0.1)
                flag_excl=1;
                code_excl=1;
            end
        end
    end
end

% A stars exclusion region

if(umag-gmag > 0.7 & umag-gmag < 1.4)     
    if(gmag-rmag > -0.5 & gmag-rmag < 0.0)
        if(rmag-imag > -0.5 & rmag-imag < 0.2) 
            if(imag-zmag > -0.4 & imag-zmag < 0.2)
                flag_excl=1;
                code_excl=2;
            end
        end
    end
end

% M stars and white dwarfs pairs exclusion region
 
if(gmag-rmag > -0.3 & gmag-rmag < 1.25)   
    if(rmag-imag > 0.6 & rmag-imag < 2.0)
        if(imag-zmag > 0.4 & imag-zmag < 1.2)
            if(errg < 0.2)
                flag_excl=1;
                code_excl=3;
            end
        end
    end
end

