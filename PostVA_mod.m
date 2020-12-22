
%This programs makes life easy by doing simple calculation on the VA
%analysed data over multiple files.
% It has a nuclear (what ever you give it to analyze) localization density
% threshold (lower and upper bound), cluster threshold (min loc. per
% cluster) and NND upper and lower distance threshhold. The filtered NND is then written in the variable 'new NND' 
% The output 'data' is used for the 'CumulativeProbDist'
% The out put file has a structyre 'VoronoiAreas', that can be used for the KL
% Divergence test written by Peter.
%

% Outputs:
% B2: A matrix with all object informations pulled from all the files
% B2: Col 1 = Cluster.nLocs
%     Col 2 = cluster.areas*pixel_to_nanometers*pixel_to_nanometers
%     Col 3 = cluster.NND*pixel_to_nanometers
%     Col 4 = Density per cluster
%     Col 5 = File index
%     Col 6 = Nuclear area corresponding to each cluster
%     Col 7 =  Cluster Areas normalized by the nuclear Area
% Obj_per_nlucleus = Matrix giving total number of Objects per nucleus for
% all the files
% nuclear_area
% obj_area_occup_perc = Matrix giving total area occupied by all the
% objects normalized by nuclear area and multiplied
% by 100, for all the files
% tot_nuclear_loc
% file_no



prompt = {'What is the pixel size of your experiment?','Enter the lower nuclear loc. density cutoff (0 if not required):','Enter higher nuclear loc. density cutoff (a very high value if not required):','Enter the Minimum Number of Localizations per Cluster (0 if not needed):', 'Enter the Maximum Number of Localizations per Cluster (a very high value if not needed):','Enter lower NND cutoff (0 if not required):','Enter higher NND cutoff (a very high value if not required):','Enter the nucleus radius cutoff in pixels:'};
title = 'Input';
dims = [1 35];
definput = {'116.99999','0', '50000000', '7', '500', '0','1000000', '300' };
answer = inputdlg(prompt,title,dims,definput);
N = str2num(answer{1});
nuclear_density_1 = str2num(answer{2});
nuclear_density_2 = str2num(answer{3});
pixel_to_nanometers = str2num(answer{1});
nloc_clust1 = str2num(answer{4});
nloc_clust2 = str2num(answer{5});
NND1=str2num(answer{6});
NND2=str2num(answer{7});
r2=str2num(answer{8});
k=0;
[fn path]=uigetfile('*.mat','Select the INPUT DATA FILE(s)','MultiSelect','on');
if iscell(fn)
    nfiles=size(fn,2);
else
    nfiles = 1;
end
clear B2;clear cargo_loc; clear tot_nuclear_loc;  clear NND_factor; clear B4; clear Obj_per_nlucleus; clear nuclear_area; clear obj_area_occup_perc;
file_no = {}; k=0; clear density_loc_per_clust; clear Obj_per_nlucleus_norm; clear tot_nuclear_loc_norm; VoronoiAreas = {}; clear data; new_NND = [];

for i=1:nfiles
    
    clear B1;clear B3; clear clust_per_cargo; clear X; clear Y; clear cargo_radius; clear Areas; clear ValidAreas;
    if nfiles == 1
        load(fullfile(path,fn));
    else
        load(fullfile(path,fn{i}));
    end
    
    X = xy(:,1); % x localizations
    Y = xy(:,2); % y localizations
    Areas = xy(:,3);
    
    %Get rid of the boundary Voronoi Areas
    K = boundary(X,Y);
    Keepers = true(length(Areas),1);
    Keepers(K) = false;
    ValidAreas = Areas(Keepers);
    
    % Calculate approximate nuclear area    
    nuclear_area_t = sum(ValidAreas)*pixel_to_nanometers*pixel_to_nanometers; % in sq nm  
    
    r = sqrt(nuclear_area_t/3.14);    
      
    % Calculate total number of localizations per cargo from raw localization
    % number
    tot_nuclear_loc_t = size(X,1);
    
        % Total nuclear loalizations normalized by nuclear area
          tot_nuclear_loc_norm_t = tot_nuclear_loc_t/nuclear_area_t;
    
  if  (tot_nuclear_loc_norm_t >= nuclear_density_1) && (tot_nuclear_loc_norm_t <= nuclear_density_2) && r<=r2  
      
    k=k+1; 
    
     if iscell(cluster)        
        B1(:,1) = cluster{1}.nLocs;
        B1(:,2) = cluster{1}.areas*pixel_to_nanometers*pixel_to_nanometers;
        B1(:,3) = cluster{1}.NND(:,1)*pixel_to_nanometers;
    else
        B1(:,1) = cluster.nLocs;
        B1(:,2) = cluster.areas*pixel_to_nanometers*pixel_to_nanometers;
        B1(:,3) = cluster.NND(:,1)*pixel_to_nanometers;
     end
    
    B1(:,4) = 0; %Density per cluster
    B1(:,5) = k;
    B1(:,6) = 0; % Crago area corresponding to each cluster
    B1(:,7) = 0; 
    
    % Removes clusters less than the threshold
    
    indices = find(B1(:,1) < nloc_clust1 | B1(:,1)> nloc_clust2 );    
    B1(indices,:) = [];
       
        
% This data is used for cumulative probability distribution

    data{k,1} = X;
    data{k,2} = Y;
    data{k,3} = Areas;
        
    % Calculate approximate cluster area    
    nuclear_area(k,1) = sum(ValidAreas)*pixel_to_nanometers*pixel_to_nanometers; % in sq nm
    
    % Calculate empty spaces per cargo, a measure of how diffused or clustered
    % the protein is. More empty spaces and higher number of clusters should represent clustered distribution.
    empty_space(k,1) =  nuclear_area(k,1)-(sum(B1(:,2)));
    
    % Ratio of cargo area to sum of cluster areas, occupied versus total area
    obj_area_occup_perc(k,1) = ((sum(B1(:,2)))/nuclear_area(k,1))*100;
      
    % Calculate total number of localizations per cargo from raw localization
    % number
    tot_nuclear_loc(1,k)= size(X,1);
    
    % Number of clusters per cargo
    s = size(B1,1);
    clust_per_cargo = s;
    
    % cargo area corresponding to each cluster
    B1(1:s,6) = nuclear_area(k,1);
    
    % Pull all the clusters from the differnt cargos  of a specific treatment
    if k==1
        B2 = B1;
        Obj_per_nlucleus = clust_per_cargo';
        VoronoiAreas{k} = ValidAreas;
    else
        B2 = vertcat(B2,B1);
        Obj_per_nlucleus = vertcat(Obj_per_nlucleus,clust_per_cargo');
        VoronoiAreas{k} = ValidAreas;
    end
    if iscell(fn)
        file_no{k} = fn{k};
    else
        file_no = fn;
    end
  end
end

if k~=0
file_no = file_no';
% Calculate density of localization per cluster
density_loc_per_clust(:,1)=B2(:,1)./B2(:,2);
B2(:,4) = density_loc_per_clust;

% Cluster Area normalized by nuclear area
B2(:,7) = B2(:,2)./B2(:,6);

tot_nuclear_loc = tot_nuclear_loc';

% Total nuclear loalizations normalized by nuclear area
tot_nuclear_loc_norm = tot_nuclear_loc./nuclear_area;

% Total number of objects per nucleus normalized by area

Obj_per_nlucleus_norm = Obj_per_nlucleus./nuclear_area;

% select NND below a threshold
indices2 = find(B2(:,3) < NND1 | B2(:,3)> NND2 );
    new_NND = B2(:,3);
    new_NND(indices2) = []; 

prompt = 'Please enter the file name: ';
filename = input(prompt,'s');
savefile = fullfile(path, [filename '.mat']);
save(savefile,'B2','nuclear_area','obj_area_occup_perc',  'file_no' , 'tot_nuclear_loc_norm', 'Obj_per_nlucleus_norm','VoronoiAreas', 'Obj_per_nlucleus', 'data', 'new_NND');
end    


