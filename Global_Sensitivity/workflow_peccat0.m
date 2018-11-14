%% Preliminaries
clear all
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';
rng(66)

% Add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 1

% a) Define the input feasible space

M = 59 ; % Number of inputs
DistrFun  = 'logn'  ; % Distribution (same for all inputs)

p0 = [0.83508546937583061; 0.69712989447674045; 10.44134948091428; ...
    9.5606753156695916; 22.263272739194772; 101.05649800089539; ...
    17.332527782983444; 19.403601740448316; 157.4047346137381; ...
    14.635262939375314; 293.37586460774401; 4.4827718526960405; ...
    300.32775190742791; 133.01397307162469; 460.66608175874353; ...
    1.1764182075409095; 84.018510931548263; 1.6867134667988213; ...
    82.485979954924986; 1.0206010478563845; 80.563272759988152; ...
    52.691607964582587;  627.25267188331998; 211.75335626665685; ...
    269.72282011062504; 1365.3428397400119; 679.52779128068778; ...
    0.62511809654037276; 1.6597626870213458; 8.5230065910935746; ...
    4.4289969350604865; 4.0997636319251773; 8.9616700444863085e-05; ...
    74.72659308067584; 0.0041468380170434221; 4.9571890876691338; ...
    3.9330111465866007; 59.735660653739174; 4.1529637851728696; ...
    9.0954018661662914; 4.4999080111784631; 1.2318431863238426; ...
    4971.2425416027772; 142.62484459355107; 0.34698445096730562; ...
    0.83843382293247926; 0.70017798153843991; 0.99091754986847957; ...
    0.51356053487506781; 0.16307907027509522; 0.50194994732634446; ...
    0.5000046943687857; 0.20186588540881914; 0.90000001580087519; ...
    0.11594685453592514; 0.90663789112969251; 0.054480065207989668; ...
    0.21671944643538985; 0.11601179754608776];

% Parameters (mu,sigma) of the lognormal distribution functions
DistrPar=cell(M,1); 
for i=1:M-15; DistrPar{i}=[log(p0(i)) abs(log(sqrt(1.5)))]; end
for i=M-14:M; DistrPar{i}=[log(p0(i)) abs(log(sqrt(p0(i))))]; end

% b) Perform sampling

%All-At-the-Time sampling (again, see workflow about RSA for more
%options):
SampStrategy = 'lhs' ;
N = 10000 ; % Number of samples
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);
[ XA, XB, XC ] = vbsa_resampling(X) ;

% c) Save to file:

% Choose name for the input file:
file_name = 'XAtight.txt' ;
file_name1 = 'XBtight.txt' ;
file_name2 = 'XCtight.txt' ;
file_name3 = 'Input_samplestight.txt' ;
% % 
% % 
% % % Set current directory to the directory where the input file will be saved:
cd([ my_dir '/example/peccat0'])
% % 
% % % Write to file:
save(file_name,'XA','-ascii') ; % With this command, data will be saved in exponential notation with 8 digits.
save(file_name1,'XB','-ascii') ;
save(file_name2,'XC','-ascii') ;
save(file_name3,'X','-ascii') ;



