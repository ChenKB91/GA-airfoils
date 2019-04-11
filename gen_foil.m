rng(42)
pop = 50;
remain = 20;
mut_rate = 0.1;
min = [0.024793 0.26 0.07 -0.7244 0.26 -0.08 0.4244 0 0 2 -15];
max = [0.024793 0.34 0.08 -0.4244 0.34 -0.07 0.7244 0 0 18 15];
%% path
par_path = [pwd '/data/para50.txt'];
force_path = [pwd '/data/force.txt'];
log_path = [pwd '/best.log'];

%% GA
% read par & force
p_fid = fopen(par_path, 'r');
f_fid = fopen(force_path,'r');
p_format = '%f %f %f %f %f %f %f %f %f %f %f\n';
f_format = '%f %f';
p_data = fscanf(p_fid,p_format,[11 50]);
f_data = fscanf(f_fid,f_format, [2 50]);

p_data = p_data.';
f_data = f_data.';
f_fitness = f_data(:,2)./f_data(:,1);
%%
% get best foils index
[~,index]  = maxk(f_fitness, 20);
best_index = index(1);
best_par   = p_data(best_index,:);
log_fid = fopen(log_path,'a');
fprintf(log_fid,'%f %f %f %f %f %f %f %f %f %f %f\n',best_par);
%%
% make pool
maxfitness = f_data(index);
prb = round(maxfitness / sum(maxfitness), 3);
%%
% select & crossover & mutation
newgen = zeros(50,11);
for i = 1:50
    ia = randsample(index,1,true,prb);
    ib = randsample(index,1,true,prb);
    par_a = p_data(ia,:);
    par_b = p_data(ib,:);
    tmp = rand(1,11)<0.5;
    child = par_a.*tmp+par_b.*~tmp;
    r_par = min + (max-min).*rand(1,11);
    r_mut = rand(1,11)<mut_rate;
    child = r_par.*r_mut + child.*~r_mut;
    newgen(i,:) = child;
end
%%
% output new param set to par_path
fclose(p_fid);
fclose(f_fid);
%%
p_fid = fopen(par_path, 'w');
fprintf(p_fid,'%f %f %f %f %f %f %f %f %f %f %f\n',newgen.');
fclose(p_fid);
%% make ib

for i = 1:50
    disp(i)
    %lastwarn('')
    
    par = newgen(i,:);
    [pts, self_cross] = evenpar(par);
    polyin = polyshape({pts(:,1 )},{pts(:,2)});
    [comx,comy] = centroid(polyin);
    s = num2str(i,' %03d');
    foil_path = [pwd '/data/ib' s '.inp'];
    [ib_fid msg] = fopen(foil_path, 'w');
    n = size(pts,1);
    % output to ib[i].inp
    % format: 1st line: npt; rest: x y coordinate
    
    if self_cross
        % write 0
        fprintf(ib_fid,'%d',0);
    else
        fprintf(ib_fid,'%d\n',n);
        fprintf(ib_fid,'%8.7f %8.7f\n',pts.');
    end
    fclose(ib_fid); 
end
exit;