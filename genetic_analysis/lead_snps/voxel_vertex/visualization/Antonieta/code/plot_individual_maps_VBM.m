
function plot_individual_maps_AP(base_data_path, SUBJECTS_DIR, base_code_dir, base_out_dir, model_name,  file_name, my_template,  my_subject, my_hemi, min_threshold, max_threshold, surf_overlay, my_surf)
base_data_path = '/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/data'; 
SUBJECTS_DIR='/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/data/SUBJECTS_DIR';
base_code_dir = '/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/code';
base_out_dir = '/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/out_plots';
model_name = 'rs41298373';%'rs7420166';%;
file_name = 'z';
my_template = 'fsaverage_sym';
min_threshold='pial';
my_subject='fake_subject';
my_hemi = 'lh'; 
min_threshold = '2';  
max_threshold = '7';  
surf_overlay = 'pial';
my_surf = 'pial';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% EXTRACTING ENVIROMENT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FREESURFER_HOME = getenv('FREESURFER_HOME');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% ADDING NECESSARY PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(fullfile(FREESURFER_HOME, 'matlab')));
addpath(genpath(base_code_dir));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% other settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
page_size = [0.25 2.5 12 9];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% LOADING FREESURFER TEMPLATE SURFACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template_surf_array = {'pial', 'inflated', 'white'};
template_full_path = fullfile (FREESURFER_HOME, 'subjects', my_template);
if exist(template_full_path, 'dir') == 7
     % * Loading template surfaces
     pial_surf = SurfStatReadSurf( {...,
            fullfile(template_full_path,'surf', 'lh.pial'),...
            fullfile(template_full_path,'surf', 'rh.pial')} );
    inflated_surf = SurfStatReadSurf( {...,
            fullfile(template_full_path,'surf', 'lh.inflated'),...
            fullfile(template_full_path,'surf', 'rh.inflated')} );   
    white_surf = SurfStatReadSurf( {...,
            fullfile(template_full_path,'surf', 'lh.white'),...
            fullfile(template_full_path,'surf', 'rh.white')} );        
             
else 
    error('template surface was not found within the FREESURFER_HOME')
    return;
end

if strcmp(surf_overlay ,'pial')
	surf=pial_surf;
elseif strcmp(surf_overlay ,'inflated')
	surf=inflated_surf;
elseif strcmp(surf_overlay ,'white')
	surf=white_surf;
end


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% CONSTRUCTING FULL PATH TO INPUT AND OUTPUT FILENAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
file_name_stripped = strrep(file_name, '.mgh' , '');
file_name_stripped = strrep(file_name_stripped, '.nii.gz' , '');
my_input_surf_texture_fname = fullfile(SUBJECTS_DIR, my_subject, 'surf', [model_name, '_', file_name_stripped , '_', my_surf,'.mgh']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% CREATING OUTPUT FOLDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(base_out_dir);
mkdir(fullfile(base_out_dir, model_name));
mkdir(fullfile(base_out_dir, model_name,  'snapshots'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% LOADING COLORMAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%load ('my_blue_white_red.mat')
load my_spectral;
load ('cmap.mat');
%load bluewhitered; 


isp('4')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% PRODUCING PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making sure the folder exists
if exist(my_input_surf_texture_fname , 'file') == 2
		    		
    %% loading texture files       
    my_texture =  load_mgh(my_input_surf_texture_fname);


    %% preparing loaded texture to be visualized with surfstat
    my_texture_tab  = transpose(squeeze(my_texture));
    if or (strcmp(my_hemi ,'lh'), strcmp(my_hemi ,'MEAN') )
        my_map_2textures = [ my_texture_tab,  zeros(size(my_texture_tab))];
    else
        my_map_2textures = [ zeros(size(my_texture_tab)),my_texture_tab];                              
    end


    %% Producing snapshots of the mean texture (area, thickness, volume, curv or sulc)
    
    my_title =  [ my_hemi,' ',  my_surf, ' ', model_name]; 
    
    my_output_fname = fullfile(base_out_dir, model_name,  'snapshots', ...
    [model_name, '_', file_name_stripped , '_', my_surf, '_', min_threshold]);
    
    
    if ~strcmp(min_threshold,'')
        my_map_2textures(abs(my_map_2textures)< abs(str2num(min_threshold)))=0;
    end
    
    
    
    [ a, cb ] = mySurfStatViewData(my_map_2textures, surf,  my_title);  
    if ~strcmp(max_threshold,'')
       SurfStatColLim( [-abs(str2num(max_threshold)), abs(str2num(max_threshold))] );
    else
        th = max(abs(my_map_2textures(:)));
        if and(th ~= 0, ~isinf(th))
            SurfStatColLim( [-th.*0.95, th.*0.95] );
        elseif (th == 0) 
            SurfStatColLim( [-1, 1] );
        else 
             th = prctile( max(abs(my_map_2textures(:))),90); 
             if and(th ~= 0, ~isinf(th))
                SurfStatColLim( [-th, th] );
             end
        end        
    end
    
    set(gca, 'FontSize',8)
     
    colormap(cmap);
    print(gcf, [my_output_fname, '_bwrCMAP.png'],'-dpng', '-r600' );
    colormap(my_spectral);
    print(gcf, [my_output_fname, '_speCMAP.png'],'-dpng', '-r600' );
end
