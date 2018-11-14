%% Compute PECCAD EET Sensitivity Indices
clear all

% Set current directory to the directory where the output file was saved:
my_dir = pwd;
addpath(genpath(my_dir))
cd([ my_dir '/example/peccat0'])

SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';    
DistrFun  = 'unif'  ; % Distribution (same for all inputs)

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

r = 10000 ; % Number of Elementary Effects 
M = 59;

X_labels = {'aMaxB', 'aMaxBP', 'aMaxF', 'KaBHiq', 'KaBLoq', 'KaBPHiq', ...
    'KaBPLoq', 'KaBPP', 'KaFHiq', 'KaFLoq', 'kBHiq', 'kBLoq', ...
    'kBPHiq', 'kBPLoq', 'kBPP', 'KdHiq', 'KdLoq', 'kFHiq', 'kFLoq', ...
    'kFP', 'KIB', 'KIF', 'kmBHiq', 'kmBLoq', 'kmBPHiq', 'kmBPLoq', ...
    'kmBPP', 'krBHiq', 'krBLoq', 'krBPHiq', 'krBPLoq', 'krBPP', ...
    'krFHiq', 'krFLoq', 'KsFP', 'mMaxB', 'mMaxBP', 'muMaxB', 'muMaxBP', ...
    'muMaxF', 'qMaxB', 'qMaxF', 'TFP', 'TyF', 'YLHiq', 'YLLoq', 'YrB', ...
    'YrF', 'YRFP', 'YsBHiq', 'YsBLoq', 'YsBPHiq', 'YsBPLoq', 'YsBPP', 'YsFHiq', ...
    'YsFLoq', 'rB0', 'rBP0', 'rF0'};

X_labels1 = {'aMaxB', 'aMaxBP', 'aMaxF', 'KaBHiq', 'KaBLoq', 'KaBPHiq', ...
    'KaBPLoq', 'KaBPP', 'KaFHiq', 'KaFLoq', 'kBHiq', 'kBLoq', ...
    'kBPHiq', 'kBPLoq', 'kBPP', 'KdHiq', 'KdLoq', 'kFHiq', 'kFLoq', ...
    'kFP', 'KIB', 'KIF', 'kmBHiq', 'kmBLoq', 'kmBPHiq', 'kmBPLoq', ...
    'kmBPP', 'krBHiq', 'krBLoq', 'krBPHiq', 'krBPLoq', 'krBPP', ...
    'krFHiq', 'krFLoq', 'KsFP', 'mMaxB', 'mMaxBP', 'muMaxB', 'muMaxBP', ...
    'muMaxF', 'qMaxB', 'qMaxF', 'TFP', 'TyF'};

xmin = p0/2;
xmax = p0*2;

DistrPar = cell(M,1);                                                          
for i=1:M; DistrPar{i} = [ xmin(i) xmax(i) ] ; end 

%% Load input and output data from files:
X = load('XMorris.txt','-ascii') ;
Y = load('YMorris.txt','-ascii') ;

%% Step 5 (Computation of the Elementary effects)                              

% Compute Elementary Effects:                                                  
[ mi, sigma ] = EET_indices(r,xmin',xmax',X,Y,design_type);     


%% Int1:  Put new and old results together                                             
r2 = 20000 ; % increase of base sample size 
X2 = load('X2Morris.txt','-ascii') ;
Ynew = load('YnewMorris.txt','-ascii') ;

Y2=[Y;Ynew]; % size (r2*(M+1),1)
   
% Compute Elementary Effects:                                                  
[ mi, sigma ] = EET_indices(r2,xmin',xmax',X2,Y2,design_type);     


%% Step 6 (Adding up new samples)                                                                                                                             
r3 = 25000 ; % increase of base sample size  
X3 = load('X3Morris.txt','-ascii') ;
Ynew1 = load('YnewMorris1.txt','-ascii') ;


Y3=[Y2;Ynew1]; % size (r2*(M+1),1)
   
% Compute Elementary Effects:                                                  
[ mi, sigma ] = EET_indices(r3,xmin',xmax',X3,Y3,design_type);

% Post-process
mi_norm = mi/max(mi);
sigma_norm = sigma/max(sigma);

% Plot L2 norms in histogram
l2norm = nan(1,size(mi_norm,2));
for i=1:size(mi_norm,2)
    l2norm(i) = norm(abs(mi_norm(1,i) - sigma_norm(1,i)));
end

% hpd results
hpdmio = [0.09, 0.05, 0.02, 0.65, 0.24, 0.26, 0.12, 0.28, 0.08, 0.06, 0.7, ...
    0.06, 0.15, 0.08, 0.77, 0.03, 0.28, 0.16, 0.02, 0.46, 0.28, 1.0, 0.13, ...
    0.06, 0.41, 0.28, 0.16, 0.13, 0.06, 0.71, 0.06, 0.38, 0.08, 0.28, 0.15, ...
    0.37, 0.34, 0.42, 0.07, 0.01, 0.10, 0.05, 0.21, 0.18, 0.67, 0.88, 1.0, ...
    0.84, 1.0, 0.23, 0.79, 0.89, 0.77, 0.73, 0.62, 0.72, 0.16, 0.29, 0.35];

% Figure
figure; 
normedl2 = l2norm/max(l2norm);
% Sort the Data & Rearrange Labels
[sorted_normedl2, new_l2indices] = sort(normedl2(1:44)); % sorts in *ascending* order
sorted_X_labels1 = X_labels1(new_l2indices);
hpd44 = hpdmio(1:44);
sorted_hpd44 = hpd44(new_l2indices);
combined = [sorted_normedl2(1:44)' ,sorted_hpd44'];
xdata = [1:44];
hb = bar(xdata, combined, 'BarWidth', 1.0);
set(gca, 'XTick', 1:length(sorted_X_labels1),'XTickLabel',sorted_X_labels1);
ylabel('Sensitivity')
%boxplot3(Si.Si(1:44))
%boxplot1(STi.STi(1:44))
xtickangle(90)
ylim([-0.05,1.15])
fn = 'Arial' ; 
fs = 8 ; 
set(gca,'FontSize',fs,'FontName',fn)
set(gca,'box','off', 'XGrid','On', 'YGrid', 'On')
set(hb(1), 'FaceColor', [153 153 153]/255)
set(hb(2), 'FaceColor', [77 77 77]/255)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.487*2, 3.487*2/1.618], 'PaperUnits', 'Inches', 'PaperSize', [3.487*2, 3.487*2/1.618])
legend_str = {'Morris Method', 'Highest Posterior Density'};
columnlegend(2,legend_str,'location','North', 'boxon');
print('GSI3','-dsvg')

