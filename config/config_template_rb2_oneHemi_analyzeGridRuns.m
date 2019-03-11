%% INITALIZE PROJECT
% evevtually it would be nice to make this into a main function to be
% called by a bash script

PROJECT_DIR = '/home/jfaskowi/JOSHSTUFF/projects/ratbrain2/' ;
mkdir(PROJECT_DIR)
cd(PROJECT_DIR)

projectFolders = { 
    'src' 
    'data'
    'bin'
} ;

for idx=1:length(projectFolders)
    addpath(genpath(strcat(PROJECT_DIR,'/',projectFolders{idx})))
end

%% SETUP GLOBAL VARS

GRID_RUN = 'results_expBern_oneHemi' ;
GRID_OUTPUT_DIR = [ PROJECT_DIR '/data/processed/grid/' GRID_RUN '/' ] ;
OUTPUT_STR = 'gridRuns1' ;

RESAMP_DIVS = cell(4,1) ;

RESAMP_DIVS{1} = struct('range', 1:11, 'size', 750, 'resamp', false) ; 
RESAMP_DIVS{2} = struct('range', 12, 'size', 500, 'resamp', false) ;
RESAMP_DIVS{3} = struct('range', 13:14, 'size', 200, 'resamp', false) ;
RESAMP_DIVS{4} = struct('range', 15:19, 'size', 100, 'resamp', true) ;

COM_NUM_RANGE = 2:20 ; 
NUM_NODES = 77 ;

REASONABLE_COM_RANGE_IND = 2:12 ;

%% make output dir

OUTPUT_DIR = strcat(PROJECT_DIR , '/data/') ; 
mkdir(OUTPUT_DIR)
