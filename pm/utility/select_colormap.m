%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function CLT = select_colormap(clt_pma, clt_select)
% obtain CLT (PMA colormap) from the enntire color lookup table
%
% OUTPUTS:
%   CLT     - Color Lookup Table
%
% INPUTS:
%   clt_pma     - Colormap Look Up Table (LUT)
%   clt_select  - Colormap Name in PMA
%
% 2/17/2019 - Yao Xiao @ SMILE | UF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CLT = select_colormap(clt_pma,clt_select)

switch clt_select
    
    case 'ASIST'
        clt_col = clt_pma.ASIST;
        
    case 'Blue_Red'
        clt_col = clt_pma.Blue_Red;
        
    case 'Gray'
        clt_col = clt_pma.Gray;
        
    case 'Inv_Gray'
        clt_col = clt_pma.Inv_Gray;
            
    case 'AJS'
        clt_col = clt_pma.AJS;
            
    case 'AZE'
        clt_col = clt_pma.AZE;
            
    case 'GE_3_Colors'
        clt_col = clt_pma.GE_3_Colors;
            
    case 'GE_Puh_Thalium'
        clt_col = clt_pma.GE_Puh_Thalium;
            
    case 'GE_Rainbow'
        clt_col = clt_pma.GE_Rainbow;
            
    case 'GE_Inv_Rainbow'
        clt_col = clt_pma.GE_Inv_Rainbow;
            
    case 'Hitachi_Block'
        clt_col = clt_pma.Hitachi_Block;
            
    case 'Hitachi_Pallete'
        clt_col = clt_pma.Hitachi_Pallete;   
        
    case 'LAS'
        clt_col = clt_pma.LAS;
                
    case 'Philips_CBF'
        clt_col = clt_pma.Philips_CBF;
                
    case 'Philips_CBV'
        clt_col = clt_pma.Philips_CBV;
                
    case 'Philips_MTT'
        clt_col = clt_pma.Philips_MTT;
                
    case 'Philips_TTP'
        clt_col = clt_pma.Philips_TTP;
                
    case 'Siemens_CT'
        clt_col = clt_pma.Siemens_CT;
                
    case 'Siemens_MR'
        clt_col = clt_pma.Siemens_MR;
                    
    case 'Terarecon'
        clt_col = clt_pma.Terarecon;
                    
    case 'Toshiba_CT_Rainbow'
        clt_col = clt_pma.Toshiba_CT_Rainbow;
                    
    case 'Toshiba_CT_Rainbow_Red'
        clt_col = clt_pma.Toshiba_CT_Rainbow_Red;
     
    case 'Toshiba_MRI'
        clt_col = clt_pma.Toshiba_MRI;
         
    case 'ZIO_NIH'
        clt_col = clt_pma.ZIO_NIH;
         
    case 'ZIO_TR'
        clt_col = clt_pma.ZIO_TR;
         
    case 'ZIO_TRW'
        clt_col = clt_pma.ZIO_TRW;
         
    case 'Tmax'
        clt_col = clt_pma.Tmax;
end

% extract CLT from the entire color lookup table
CLT = zeros(256,3);
CLT(:,1) = clt_col(1:256);
CLT(:,2) = clt_col(259:514);
CLT(:,3) = clt_col(517:772);
CLT = CLT./255;

end

