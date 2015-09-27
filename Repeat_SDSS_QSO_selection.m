%Please cite Murphy, M.T., Bernet, M.L., 2015, MNRAS, submitted, if you use this code
%Comments or reports of errors are welcome. Email: bernet.martin@gmail.com
%
%Please check the readme.txt file for the applicability of the algorithm and its
%restrictions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we provide the data (List_QSO_paperMMMB_TARGETmag_coord.txt), which we used in 
%(Murphy, M.T., Bernet, M.L., 2015, MNRAS, submitted) to test for a dust bias

%The file List_QSO_DR7_total_TARGET.txt contains the full SDSS QSO DR7 sample
%as published in Schneider D. P. et al., 2010, AJ, 139, 2360



[Au,mjd,plate,fiber,zQSO_t,Lzf,Hzf,Firstf,Rosatf,Serendipitysf,Starf,Galf,umag_t,erru_t,gmag_t,errg_t,rmag_t,errr_t,imag_t,erri_t,zmag_t,errz_t,ra_t,dec_t,sdss_t]=textread('./List_QSO_paperMMMB_TARGETmag_coord.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s'); 

%[Au,mjd,plate,fiber,zQSO_t,Lzf,Hzf,Firstf,Rosatf,Serendipitysf,Starf,Galf,umag_t,erru_t,gmag_t,errg_t,rmag_t,errr_t,imag_t,erri_t,zmag_t,errz_t,ra_t,dec_t,sdss_t]=textread('./List_QSO_DR7_total_TARGET.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s'); 


dr7_sch(:,1)=umag_t;
dr7_sch(:,2)=gmag_t;
dr7_sch(:,3)=rmag_t;
dr7_sch(:,4)=imag_t;
dr7_sch(:,5)=zmag_t;
dr7_sch(:,6)=erru_t;
dr7_sch(:,7)=errg_t;
dr7_sch(:,8)=errr_t;
dr7_sch(:,9)=erri_t;
dr7_sch(:,10)=errz_t;
dr7_sch(:,11)=Au;         %Au 
EBmV=Au./5.155;           % E(B-V) Galactic extinction      
dr7_sch(:,12)=3.793*EBmV; %Ag
dr7_sch(:,13)=2.751*EBmV; %Ar
dr7_sch(:,14)=2.086*EBmV; %Ai
dr7_sch(:,15)=1.479*EBmV; %Az

%Correct for Galactic extinction

dr7_sch(:,1)=dr7_sch(:,1)-dr7_sch(:,11);  
dr7_sch(:,2)=dr7_sch(:,2)-dr7_sch(:,12);
dr7_sch(:,3)=dr7_sch(:,3)-dr7_sch(:,13);
dr7_sch(:,4)=dr7_sch(:,4)-dr7_sch(:,14);
dr7_sch(:,5)=dr7_sch(:,5)-dr7_sch(:,15);
        
redvector_u(:,1)=10.0.^(0.4.*dr7_sch(:,11));
redvector_g(:,1)=10.0.^(0.4.*dr7_sch(:,12));
redvector_r(:,1)=10.0.^(0.4.*dr7_sch(:,13));
redvector_i(:,1)=10.0.^(0.4.*dr7_sch(:,14));
redvector_z(:,1)=10.0.^(0.4.*dr7_sch(:,15));
        
dr7_sch(:,6)=dr7_sch(:,6)./sqrt(redvector_u(:,1));
dr7_sch(:,7)=dr7_sch(:,7)./sqrt(redvector_g(:,1));
dr7_sch(:,8)=dr7_sch(:,8)./sqrt(redvector_r(:,1));
dr7_sch(:,9)=dr7_sch(:,9)./sqrt(redvector_i(:,1));
dr7_sch(:,10)=dr7_sch(:,10)./sqrt(redvector_z(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
dr7_sch(:,17)=zQSO_t;
dr7_sch(:,18)=ra_t;
dr7_sch(:,19)=dec_t;

%Flags according to Schneider et al. 2010 (QSO catalogue DR7)

dr7_sch(:,20)=Lzf;
dr7_sch(:,21)=Hzf;
dr7_sch(:,22)=Firstf;
dr7_sch(:,23)=Rosatf;
dr7_sch(:,24)=Serendipitysf;
dr7_sch(:,25)=Starf;
dr7_sch(:,26)=Galf;
dr7_sch(:,28)=EBmV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do ugri color selection   

i_lowz_outl=1;
i_lowz_incl=1;
i_color_comb=1;
i_highz_outl=1;
i_highz_incl=1;
i_incl_comb=1;
i_outl_comb=1;

i_lowz_schn=1;
i_highz_schn=1;
i_firstf_schn=1;
i_rosatf_schn=1;
i_serend_schn=1;
i_galf_schn=1;
i_starf_schn=1;
i_comb_schn=1;
i_midz=1;
i_total=1;


for i=1:length(dr7_sch(:,1))
    
     [dr7_sch(i,29),temp1,temp2,dr7_sch(i,33)]=stellar_locus_ugri(dr7_sch(i,1),dr7_sch(i,2),dr7_sch(i,3),dr7_sch(i,4),dr7_sch(i,6),dr7_sch(i,7),dr7_sch(i,8),dr7_sch(i,9),4.0,19.1);
     
     [dr7_sch(i,30),temp3,temp4,dr7_sch(i,34)]=stellar_locus_griz(dr7_sch(i,1),dr7_sch(i,2),dr7_sch(i,3),dr7_sch(i,4),dr7_sch(i,5),dr7_sch(i,7),dr7_sch(i,8),dr7_sch(i,9),dr7_sch(i,10),4.0);
     %flag low_z    high_z
     [dr7_sch(i,31),dr7_sch(i,32),dr7_sch(i,35)]=flagging_incl_region(dr7_sch(i,1),dr7_sch(i,2),dr7_sch(i,3),dr7_sch(i,4),dr7_sch(i,5),dr7_sch(i,6),dr7_sch(i,7),dr7_sch(i,8),dr7_sch(i,9));
     
     [dr7_sch(i,36),dr7_sch(i,37)]=flagging_excl_region(dr7_sch(i,1),dr7_sch(i,2),dr7_sch(i,3),dr7_sch(i,4),dr7_sch(i,5),dr7_sch(i,6),dr7_sch(i,7),dr7_sch(i,8),dr7_sch(i,9),dr7_sch(i,10));
     
     if((dr7_sch(i,29)==1 | dr7_sch(i,30)==1 | dr7_sch(i,31)==1 | dr7_sch(i,32)==1) & dr7_sch(i,36)==0)
         
         dr7_total_comb(i_total,:)=dr7_sch(i,:);
         i_total=i_total+1;
     end  
     
     if(dr7_sch(i,29)==1 & dr7_sch(i,36)==0)
         
         dr7_lowz_outl(i_lowz_outl,:)=dr7_sch(i,:);
         i_lowz_outl=i_lowz_outl+1;
     end
     
     if(dr7_sch(i,31)==1 & dr7_sch(i,36)==0)
         
         dr7_lowz_incl(i_lowz_incl,:)=dr7_sch(i,:);
         i_lowz_incl=i_lowz_incl+1;
     end
     
     if((dr7_sch(i,29)==1 | dr7_sch(i,30)==1) & dr7_sch(i,36)==0)
         
         dr7_color_comb(i_color_comb,:)=dr7_sch(i,:);
         i_color_comb=i_color_comb+1;
     end
     
     if((dr7_sch(i,31)==1 | dr7_sch(i,32)==1) & dr7_sch(i,36)==0)
         
         dr7_incl_comb(i_incl_comb,:)=dr7_sch(i,:);
         i_incl_comb=i_incl_comb+1;
     end
     
     if(dr7_sch(i,30)==1 & dr7_sch(i,36)==0)
         
         dr7_highz_outl(i_highz_outl,:)=dr7_sch(i,:);
         i_highz_outl=i_highz_outl+1;
     end
     
     if(dr7_sch(i,32)==1 & dr7_sch(i,36)==0)
         
         dr7_highz_incl(i_highz_incl,:)=dr7_sch(i,:);
         i_highz_incl=i_highz_incl+1;
     end
     
     if(dr7_sch(i,20)==1)
         
         dr7_lowz_schn(i_lowz_schn,:)=dr7_sch(i,:);
         i_lowz_schn=i_lowz_schn+1;
     end
     
     if(dr7_sch(i,21)==1)
         
         dr7_highz_schn(i_highz_schn,:)=dr7_sch(i,:);
         i_highz_schn=i_highz_schn+1;
     end
     
     if(dr7_sch(i,20)==1 | dr7_sch(i,21)==1)
         
         dr7_comb_schn(i_comb_schn,:)=dr7_sch(i,:);
         i_comb_schn=i_comb_schn+1;
     end
     
     if(dr7_sch(i,24)==1)
         
         dr7_serend_schn(i_serend_schn,:)=dr7_sch(i,:);
         i_serend_schn=i_serend_schn+1;
     end
         
end



figure(1)
subplot(2,1,1)
plot(dr7_lowz_schn(:,1)-dr7_lowz_schn(:,2),dr7_lowz_schn(:,2)-dr7_lowz_schn(:,3),'ok')
xlabel('u-g','fontsize',18,'fontweight','b');
ylabel('g-r','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_lowz_schn(:,2)-dr7_lowz_schn(:,3),dr7_lowz_schn(:,3)-dr7_lowz_schn(:,4),'ok')
xlabel('g-r','fontsize',18,'fontweight','b')
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

figure(2)
subplot(2,1,1)
plot(dr7_lowz_outl(:,1)-dr7_lowz_outl(:,2),dr7_lowz_outl(:,2)-dr7_lowz_outl(:,3),'ob')
hold on
plot(dr7_lowz_incl(:,1)-dr7_lowz_incl(:,2),dr7_lowz_incl(:,2)-dr7_lowz_incl(:,3),'or')

xlabel('u-g','fontsize',18,'fontweight','b');
ylabel('g-r','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_lowz_outl(:,2)-dr7_lowz_outl(:,3),dr7_lowz_outl(:,3)-dr7_lowz_outl(:,4),'ob')
hold on
plot(dr7_lowz_incl(:,2)-dr7_lowz_incl(:,3),dr7_lowz_incl(:,3)-dr7_lowz_incl(:,4),'or')
xlabel('g-r','fontsize',18,'fontweight','b')
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);


figure(3)
subplot(2,1,1)
plot(dr7_highz_schn(:,2)-dr7_highz_schn(:,3),dr7_highz_schn(:,3)-dr7_highz_schn(:,4),'ok')
xlabel('g-r','fontsize',18,'fontweight','b');
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_highz_schn(:,3)-dr7_highz_schn(:,4),dr7_highz_schn(:,4)-dr7_highz_schn(:,5),'ok')
xlabel('r-i','fontsize',18,'fontweight','b')
ylabel('i-z','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

figure(4)
subplot(2,1,1)
plot(dr7_highz_outl(:,2)-dr7_highz_outl(:,3),dr7_highz_outl(:,3)-dr7_highz_outl(:,4),'ob')
hold on
plot(dr7_highz_incl(:,2)-dr7_highz_incl(:,3),dr7_highz_incl(:,3)-dr7_highz_incl(:,4),'or')

xlabel('g-r','fontsize',18,'fontweight','b');
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_highz_outl(:,3)-dr7_highz_outl(:,4),dr7_highz_outl(:,4)-dr7_highz_outl(:,5),'ob')
hold on
plot(dr7_highz_incl(:,3)-dr7_highz_incl(:,4),dr7_highz_incl(:,4)-dr7_highz_incl(:,5),'or')
xlabel('r-i','fontsize',18,'fontweight','b')
ylabel('i-z','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

figure(5)

subplot(2,1,1)
plot(dr7_comb_schn(:,1)-dr7_comb_schn(:,2),dr7_comb_schn(:,2)-dr7_comb_schn(:,3),'ok')
xlabel('u-g','fontsize',18,'fontweight','b');
ylabel('g-r','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_comb_schn(:,2)-dr7_comb_schn(:,3),dr7_comb_schn(:,3)-dr7_comb_schn(:,4),'ok')
xlabel('g-r','fontsize',18,'fontweight','b')
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

figure(6)

subplot(2,1,1)
plot(dr7_total_comb(:,1)-dr7_total_comb(:,2),dr7_total_comb(:,2)-dr7_total_comb(:,3),'ok')
xlabel('u-g','fontsize',18,'fontweight','b');
ylabel('g-r','fontsize',18,'fontweight','b')
axis([-0.5,3.2,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);

subplot(2,1,2)
plot(dr7_total_comb(:,2)-dr7_total_comb(:,3),dr7_total_comb(:,3)-dr7_total_comb(:,4),'ok')
xlabel('g-r','fontsize',18,'fontweight','b')
ylabel('r-i','fontsize',18,'fontweight','b')
axis([-0.5,3.0,-0.5,3.0]);
grid(gca,'minor')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [0.025,0.06] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'FontWeight'   ,'Bold'         ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'Layer','Top'             , ...
  'FontSize',18);



