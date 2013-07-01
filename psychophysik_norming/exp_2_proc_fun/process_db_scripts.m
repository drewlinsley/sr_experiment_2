function process_db_scripts(settings)
%produces a stim X block X run matrix


if ~exist('settings','var'),
    settings.processed_dbs_dir = fullfile(out_dir,'processed_dbs');
    settings.in_dir = '/Users/drewlinsley/Documents/Dropbox/deBruijn_MacOSX';
    settings.trash_dir = fullfile(settings.in_dir,'db_output');
    settings.out_name = 'current_sub';
    settings.out_dir = fullfile(settings.in_dir,'exp_2_scripts',settings.out_name);
    settings.script_dir = fullfile(out_dir,out_name);
    settings.k = 6;
    settings.n = 3;
    settings.num_blocks = 6;
    settings.numbins = 6;
    settings.distfile = '/Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/exp_2_distances.txt';
    settings.isi = 1500;
    settings.tmin = 10;
    settings.tmax = 79;
    settings.numiter = 50;
    settings.seshlimit = 8;
    settings.pad_seq = 'off';
    settings.wrap_seq = 'off';
end
if ~exist(settings.processed_dbs_dir,'dir'),
    mkdir(settings.processed_dbs_dir)
end
num_stims = (settings.k^settings.n)+((settings.k^settings.n)/settings.k);

in_files = dir(fullfile(settings.out_dir,sprintf('%s*',settings.out_name)));

db_array = nan(num_stims,settings.seshlimit,settings.num_runs);
file_counter = 1;
for r_idx = 1:settings.num_runs,
    for b_idx = 1:settings.seshlimit,
        fid = fopen(fullfile(settings.out_dir,in_files(file_counter).name));
        file_data = textscan(fid,'%s','Delimiter',':');
        fclose(fid);
        if isnan(str2double(file_data{1}(3,:))),
            file_data = file_data{1}(6,:);
        else
            file_data = file_data{1}(3,:);
        end
        file_data = strrep(file_data,'0','00'); %all nulls are double time
        db_array(:,b_idx,r_idx) = arrayfun(@str2double,file_data{1})';
        file_counter = file_counter + 1;
    end
end

switch settings.wrap_seq
    case 'on'
        end_seqs = db_array(numel(db_array(:,1))-9:numel(db_array(:,1)),numel(db_array(1,:))); %take off last 10 entries of last run
        end_seqs = [end_seqs,db_array(numel(db_array(:,1))-9:numel(db_array(:,1)),1:numel(db_array(1,:))-1)]; %add to the front of matrix of the last 10 entries of other runs
        db_array = [end_seqs;db_array]; %put the sandwich together, with 10 null events at the end
end

switch settings.pad_seq
    case 'on'
        db_array = [db_array;zeros(10,settings.seshlimit,settings.num_runs)]; %put the sandwich together, with 10 null events at the end
end


save(fullfile(settings.processed_dbs_dir,strcat(settings.out_name,'.mat')),'db_array');