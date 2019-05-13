function gen_foil(now_gen)
rng(42)
pop = 70;
remain = 30;
mut_rate = 0.07;
% NACA0015
% 0.024793 0.300000 0.075000 -0.574473 0.300000 -0.075000 0.578907 0.000000 0.000000 9.954040 -5.254758
min = [0.024793 0.25 0.07 -0.8244 0.26 -0.08 0.4244 -0.05 0.015 8 -60];
max = [0.024793 0.35 0.08 -0.4244 0.34 -0.07 0.8244 0.05 0.015 12 60];
%% path
par_path = [pwd '/data/para50.txt'];
force_path = [pwd '/data/force.txt'];

best_path = [pwd '/log/best.log'];
mid_path = [pwd '/log/mid.log'];
worst_path = [pwd '/log/worst.log'];

best_force_path = [pwd '/log/best_force.log'];
mid_force_path = [pwd '/log/mid_force.log'];
worst_force_path = [pwd '/log/worst_force.log'];

%% GA
% read par & force
p_fid = fopen(par_path, 'r');
f_fid = fopen(force_path,'r');
p_format = '%f %f %f %f %f %f %f %f %f %f %f\n';
f_format = '%f %f';
p_data = fscanf(p_fid,p_format,[11 pop]);
f_data = fscanf(f_fid,f_format, [2 pop]);

p_data = p_data.';
f_data = f_data.';
f_fitness = f_data(:,2)./f_data(:,1);
%%
% get best foils index
[~,sorted_index]  = sort(f_fitness);
sorted_index = flipud(sorted_index);
index = sorted_index(1:20);

best_index = index(1);
mid_index = sorted_index(pop/2);
worst_index = sorted_index(pop);

best_par   = p_data(best_index,:);
mid_par   = p_data(mid_index,:);
worst_par   = p_data(worst_index,:);

log_fid = fopen(best_path,'a');
flog_fid = fopen(best_force_path,'a');
disp(best_par)
fprintf(log_fid,'%d %f %f %f %f %f %f %f %f %f %f %f\n',[now_gen best_par]);
fprintf(flog_fid,'%d, %f, %f, %f\n',[now_gen f_data(best_index,:) f_fitness(best_index)]);

log2_fid = fopen(mid_path,'a');
flog2_fid = fopen(mid_force_path,'a');
disp(mid_par)
fprintf(log2_fid,'%d %f %f %f %f %f %f %f %f %f %f %f\n',[now_gen mid_par]);
fprintf(flog2_fid,'%d, %f, %f, %f\n',[now_gen f_data(mid_index,:) f_fitness(mid_index)]);

log3_fid = fopen(worst_path,'a');
flog3_fid = fopen(worst_force_path,'a');
disp(worst_par)
fprintf(log3_fid,'%d %f %f %f %f %f %f %f %f %f %f %f\n',[now_gen worst_par]);
fprintf(flog3_fid,'%d, %f, %f, %f\n',[now_gen f_data(worst_index,:) f_fitness(worst_index)]);
%%
% make pool
maxfitness = f_fitness(index);
prb = round(maxfitness / sum(maxfitness), 3);
%%
% select & crossover & mutation
newgen = zeros(pop,11);
remain = 3;
for i = 1:remain
    newgen(i,:) = p_data(index(i),:);
end
for i = remain+1:pop
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

for i = 1:pop
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
        disp('self-cross... writing 0')
        fprintf(ib_fid,'%d',0);
    else
        fprintf(ib_fid,'%d\n',n);
        fprintf(ib_fid,'%8.7f %8.7f\n',pts.');
    end
    fclose(ib_fid);
end

disp(['Done Processing generation ' num2str(now_gen)])
% exit;