clear
close all
clc
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
myfile = 'oddball_newsubs_3.csv';
study_name = 'BC_scene_classification';%pilot_classification or pilot_classification_with_training
produce_fmri_fits = 'on';
generate_dbs_for_each = 'on';
%--- After running with the above settings on, you can turn the below
%settings on -- of course, process_dbs only works after executing their db
%scripts
dbs_to_execute = [2,4,5,6,9,10,12,13,16,17,18];
execute_dbs_for_each = 'off';
process_dbs = 'off';
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
min_rooms_per_grouping = 50; %average bathroom/kitchen/neutral groupings must have this many rooms in them
low_group = 0.333; %average bathrooms start point
mid_group = 0.5; %neutral point
high_group = 0.667; %average kitchens start point
opts = statset('MaxIter',1000); %iterations for curve fit
%-- plotting params
x_plots = 3;
y_plots = 4;
%%%%%%%%%%%%%%%%%%%
%these parameters must be changed for each study
home_dir = '/Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming';
in_dir = fullfile(home_dir,'psychophysik_datasheets');
stims_to_toss = 20; %if no training is used, toss at least 1; otherwise toss the # of training stims
group_width = 50;
total_images = 300;
toss_header = 5;
toss_footer = 6;%OS/Browser/Xres/YRes/Date/Time
num_stims = 120;
image_groupings = group_width:group_width:total_images;
s_image_groupings = 1:group_width:total_images;
%%%%%%%%%%%%%%%%%%%
settings.in_dir = '/Users/drewlinsley/Documents/Dropbox/deBruijn_MacOSX'; %location of db program
settings.k = 6; %num stims
settings.n = 2; %amount of counterbalancing
settings.numbins = 6;
settings.num_runs = 8;
settings.isi = 1500;
settings.tmin = 10;
settings.tmax = 79;
settings.numiter = 50;
settings.seshlimit = 6; %number of blocks
settings.nice = 19;% set script niceness -20 -> 20 with lower  valshaving higher priority
settings.pad_seq = 'off';
settings.wrap_seq = 'off';
settings.trash_dir = fullfile(home_dir,'db_output');
settings.script_dir = fullfile(home_dir,'exp_2_scripts');
settings.distfile = fullfile(home_dir,'exp_2_distances.txt');
settings.processed_dbs_dir = fullfile(home_dir,'processed_dbs');
%%%%%%%%%%%%%%%%%%%
%these are optional
preprocess = 'rt_log_trim';%rt_trim (adaptive) or fixed_rt
rt_modifier = 4; %for adaptive RT
rt_floor = 0.2; %for fixed RT
rt_ceiling = 5.0; %for fixed RT
feature_selection = '';
timing_trim_param = 4; %#SDs
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%l
fid = fopen(myfile);
file_data = textscan(fid,'%s','Delimiter','"');
fclose(fid);
file_data = file_data{1};
%file_data = importdata(fullfile(in_dir,myfile));

rt = nan(numel(file_data(:,1)),num_stims);
data = nan(numel(file_data(:,1)),num_stims);
key = nan(numel(file_data(:,1)),num_stims);
timing = nan(numel(file_data(:,1)),num_stims);
sub_emails = cell(numel(file_data(:,1)),1);

key_b = nan(numel(file_data(:,1)),num_stims);

key_trimmer_a = @(x) strrep(x,sprintf('/psychophysik/%s/',study_name),'');
key_trimmer_b = @(x) strrep(x,'.jpg','');
key_trimmer_c = @(x) strrep(x,'b_','');
key_trimmer_d = @(x) strrep(x,'k_','');

for sub_idx = 1:numel(data(:,1)),
    this_data = file_data{sub_idx};
    this_data = regexp(this_data,',','split');
    sub_emails{sub_idx} = this_data{4}; %save subject email
    t_this_data = this_data(:,toss_header+1:(numel(this_data))-toss_footer);
    results = t_this_data(1:num_stims);
    stim_time = t_this_data((1+num_stims):(num_stims*2));
    reaction_time = t_this_data(((num_stims*2)+1):(num_stims*3));
    stim_key = t_this_data(((num_stims*3)+2):((num_stims*4)+1));
    
    
    t_data = strrep(results,'Kitchen','1');
    t_data = strrep(t_data,'Bathroom','0');
    data(sub_idx,:) = str2double(t_data);
    timing(sub_idx,:) = str2double(stim_time);
    rt(sub_idx,:) = str2double(reaction_time);
    t_key = cellfun(key_trimmer_a,stim_key,'UniformOutput',false);
    t_key = cellfun(key_trimmer_b,t_key,'UniformOutput',false);
    t_key = cellfun(key_trimmer_c,t_key,'UniformOutput',false);
    t_key = cellfun(key_trimmer_d,t_key,'UniformOutput',false);
    key(sub_idx,:) = str2double(t_key);
end
%--- Fix email addresses
fix_email = @(x) strrep(x,'_','\_');
sub_emails = cellfun(fix_email,sub_emails,'UniformOutput',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save training runs
training_rt = rt(:,(2:(stims_to_toss)));
training_data = data(:,(2:(stims_to_toss)));
training_key = key(:,(2:(stims_to_toss)));
training_timing = timing(:,(2:(stims_to_toss)));
%process training accuracy
training_key(training_key<151)=0;
training_key(training_key>0)=1;
training_accuracy = round(sum(training_data==training_key,2)./numel(training_data(1,:))*100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process actual data
rt = rt(:,((stims_to_toss+1):(numel(rt(1,:))-1)));
data = data(:,((stims_to_toss+1):(numel(data(1,:))-1)));
key = key(:,((stims_to_toss+1):(numel(key(1,:))-1)));
timing = timing(:,((stims_to_toss+1):(numel(timing(1,:))-1)));

%remove any negative timing/rts... consider these program errors
data(timing<0)=NaN;
rt(timing<0)=NaN;
key(timing<0)=NaN;
timing(timing<0)=NaN;
data(rt<0)=NaN;
key(rt<0)=NaN;
timing(rt<0)=NaN;
rt(rt<0)=NaN;

rt = rt / 1000;
timing = timing / 1000;


%trim timing for + timing_trim_paramSDs
timing_m = mean(mean(timing,2));
timing_sd = std(reshape(timing,numel(timing),1));
data(timing>=(timing_m + (timing_trim_param*timing_sd)))=NaN;
key(timing>=(timing_m + (timing_trim_param*timing_sd)))=NaN;
rt(timing>=(timing_m + (timing_trim_param*timing_sd)))=NaN;

%%%%%%create psychometric function
switch preprocess
    case 'rt_trim'
        rt_sd = nanstd(reshape(rt,numel(rt),1));
        rt_mean = nanmean(nanmean(rt,2));
        data(rt<=(rt_mean-(rt_sd*rt_modifier))|rt>=(rt_mean+(rt_sd*rt_modifier)))=NaN;
        fprintf('\rTrimmed Data for RT excursions')
    case 'rt_log_trim'
        l_rt = log(rt);
        rt_sd = nanstd(reshape(l_rt,numel(rt),1));
        rt_mean = nanmean(nanmean(l_rt,2));
        data(l_rt<=(rt_mean-(rt_sd*rt_modifier))|l_rt>=(rt_mean+(rt_sd*rt_modifier)))=NaN;
        fprintf('\rTrimmed Data for RT excursions')
    case 'fixed_rt'
        data(rt<=(rt_floor)|rt>=(rt_ceiling))=NaN;
        fprintf('\rFix-Trimmed Data for RT excursions')
    otherwise
end
switch feature_selection
    case 'maximize'
        counter = 1;
        minimum_samples = 40;
        maximum_samples = numel(data(1,:));
        t_mat = nan(maximum_samples-minimum_samples,1);
        for curr_ind = minimum_samples:maximum_samples,
            t_data = data(:,(1:curr_ind));
            t_key = key(:,(1:curr_ind));
            t_key(t_key<=150) = 0;
            t_key(t_key~=0)=1;
            t_data = t_data + t_key; %incorrect answers are 1s
            proportions = nan(numel(t_data(:,1)),1);
            for mm = 1:numel(proportions),
                proportions(mm) = sum(t_data(mm,:)==0 | t_data(mm,:)==2)/ numel(t_data(mm,:));
            end
            t_mat(counter) = (nanmean(proportions)-.5)/(nanstd(proportions)/sqrt(numel(proportions)));
            counter = counter + 1;
        end
        t_mat(:,2) = minimum_samples:maximum_samples;
        end_point = t_mat(t_mat(:,1)==max(t_mat(:,1)),2);
        rt = rt(:,1:end_point);
        data = data(:,1:end_point);
        key = key(:,1:end_point);
        fprintf('\rChose the sample size that Maximized Difference from Chance\r')
    otherwise
end
bathroom_mat = nan(numel(image_groupings(:,1)),numel(data(:,1))); %bathroom intensity X subject matrix
kitchen_mat = nan(numel(image_groupings(:,1)),numel(data(:,1))); %kitchen intensity X subject matrix
num_judgments = nan(numel(image_groupings(:,1)),numel(data(:,1))); %kitchen intensity X subject matrix

correct_mat = nan(numel(image_groupings(:,1)),numel(data(:,1))); %kitchen intensity X subject matrix
incorrect_mat = nan(numel(image_groupings(:,1)),numel(data(:,1)));
intensity_ind = nan(numel(image_groupings(:,1)),numel(data(:,1)));
for subject = 1:numel(data(:,1)),
    for mm = 1:numel(image_groupings),
        low_ind = image_groupings(mm)-(group_width-1);
        high_ind = image_groupings(mm);
        selected_imgs = key(subject,key(subject,:)>=low_ind&key(subject,:)<=high_ind);
        selected_imgs = selected_imgs(~isnan(selected_imgs));
        selected_judgments = data(subject,key(subject,:)>=low_ind&key(subject,:)<=high_ind); %bathrooms are 0s, kitchens are 1s
        selected_judgments = selected_judgments(~isnan(selected_judgments));
        if image_groupings(mm)<=(total_images/2), %bathroom is correct
            correct_mat(mm,subject) = sum(selected_judgments==0);
            incorrect_mat(mm,subject) = sum(selected_judgments==1);
        else %kitchen is correct
            correct_mat(mm,subject) = sum(selected_judgments==1);
            incorrect_mat(mm,subject) = sum(selected_judgments==0);
        end
        bathroom_mat(mm,subject) = sum(selected_judgments==0);
        kitchen_mat(mm,subject) = sum(selected_judgments==1);
        num_judgments(mm,subject) = numel(selected_judgments);
        intensity_ind(mm,subject) = ((image_groupings(mm)-(group_width/2))-(total_images/2));
    end
end

if sum(kitchen_mat)==0,
    warning('I''m guessing you didn''t change "study name"... Annoying, I know')
end

prop_mat = correct_mat./(incorrect_mat+correct_mat);
prop_bathroom_mat = bathroom_mat ./ num_judgments;
prop_kitchen_mat = kitchen_mat ./ num_judgments;

intensity_ind = (intensity_ind-min(intensity_ind(:,1)))./(max(intensity_ind(:,1))-min(intensity_ind(:,1)));

%%%interpolate nans
%remove naned out proportion data
counter = 1;
for mm = 1:numel(prop_mat(1,:)),
    if(sum(isnan(prop_mat(:,counter)))>=(numel(prop_mat(:,1)))/2), %if at least half of a participant's runs are NaNed
        t_prop = prop_mat(:,~ismembc(1:numel(prop_mat(1,:)),counter));
        prop_mat = t_prop;
    else
        counter = counter + 1;
    end
end

%interpolate
if sum(sum(isnan(prop_mat)))~=0,
    for mm = 1:numel(prop_mat(1,:)),
        %propmat nans
        if isnan(prop_mat(1,mm)),
            prop_mat(1,mm) = prop_mat(find(~isnan(prop_mat(:,mm)), 1 ),mm);
        end
        if isnan(prop_mat(numel(prop_mat(:,1)),mm)),
            prop_mat(numel(prop_mat(:,1)),mm) = prop_mat(find(~isnan(prop_mat(:,mm)), 1, 'last' ),mm);
        end
        nanx = isnan(prop_mat(:,mm));
        t = 1:numel(prop_mat(:,mm));
        prop_mat(nanx,mm) = interpft(t(~nanx), prop_mat((~nanx),mm), t(nanx));
        %prop bathrooms nans
        if isnan(prop_bathroom_mat(1,mm)),
            prop_bathroom_mat(1,mm) = prop_bathroom_mat(find(~isnan(prop_bathroom_mat(:,mm)), 1 ),mm);
        end
        if isnan(prop_bathroom_mat(numel(prop_bathroom_mat(:,1)),mm)),
            prop_bathroom_mat(numel(prop_bathroom_mat(:,1)),mm) = prop_bathroom_mat(find(~isnan(prop_bathroom_mat(:,mm)), 1, 'last' ),mm);
        end
        nanx = isnan(prop_bathroom_mat(:,mm));
        t = 1:numel(prop_bathroom_mat(:,mm));
        prop_bathroom_mat(nanx,mm) = interpft(t(~nanx), prop_bathroom_mat((~nanx),mm), t(nanx));
        %prop kitchen nans
        if isnan(prop_kitchen_mat(1,mm)),
            prop_kitchen_mat(1,mm) = prop_kitchen_mat(find(~isnan(prop_kitchen_mat(:,mm)), 1 ),mm);
        end
        if isnan(prop_kitchen_mat(numel(prop_kitchen_mat(:,1)),mm)),
            prop_kitchen_mat(numel(prop_kitchen_mat(:,1)),mm) = prop_kitchen_mat(find(~isnan(prop_kitchen_mat(:,mm)), 1, 'last' ),mm);
        end
        nanx = isnan(prop_kitchen_mat(:,mm));
        t = 1:numel(prop_kitchen_mat(:,mm));
        prop_kitchen_mat(nanx,mm) = interpft(t(~nanx), prop_kitchen_mat((~nanx),mm), t(nanx));
    end
end


%---Toss any dudes who show 0 variability... didn't understand the task
subs_tossed = sum(std(prop_kitchen_mat)==0);
prop_kitchen_mat = prop_kitchen_mat(:,(std(prop_kitchen_mat)>0));
%---Flip the props for any dude who reversed buttons (i.e. a negative
%slope)
for sub_idx = 1:numel(prop_kitchen_mat(1,:)),
    fit = regstats(prop_kitchen_mat(:,1),(1:total_images/group_width)','linear',{'tstat','beta'});
    if fit.tstat.pval(2) < 0.05,
        prop_kitchen_mat(:,sub_idx) = prop_kitchen_mat(numel(prop_kitchen_mat(:,sub_idx)):-1:1,sub_idx);
        fprintf('\r\rFlipped %s\r',sub_emails{sub_idx})
    end
end

%%%%%%%
%proportion abs kitchens
counter = 1;
for aa = 1:numel(prop_kitchen_mat(1,:)),
    if counter == 1,
        figure,
        hold on
    end
    subplot(x_plots,y_plots,counter)
    %[curve goodness] = fit(abs(intensity_ind(:,aa)),prop_kitchen_mat(:,aa),'cubicspline');
    %plot(curve,abs(intensity_ind(:,aa)),(prop_kitchen_mat(:,aa)),'.')
    plot(intensity_ind(:,aa),(prop_kitchen_mat(:,aa)),'o')
    
    legend('off')
    xlabel('Proportion of Kitchen Spaciousness');
    ylabel('Proportion Kitchen Judgments');
    title({sprintf('S%i Training accuracy = %i%%',aa,training_accuracy(aa)),sprintf('%s',sub_emails{aa})})
    axis([-0.1 1.1 -0.1 1.1])
    counter = counter + 1;
    if counter > x_plots*y_plots,
        fig_title = sprintf('Subject Performance -- %i training trials and %i Subjects tossed',stims_to_toss,subs_tossed);
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', fig_title, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        counter = 1;
    end
end

fig_title = sprintf('Subject Performance -- %i training trials and %i Subjects tossed',stims_to_toss,subs_tossed);
annotation('textbox', [0 0.9 1 0.1], ...
    'String', fig_title, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

%--- Fit weibull curves to find neutral point
counter = 1;
my_intensity_ind = intensity_ind + 0.001; %for curve fitting
sub_mids = nan(numel(prop_kitchen_mat(1,:)),1);
sub_average_kitchens = nan(numel(prop_kitchen_mat(1,:)),1);
sub_average_bathrooms = nan(numel(prop_kitchen_mat(1,:)),1);
f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
min_prop = (min_rooms_per_grouping*3)/total_images; %minimum proportion spanned by sigmoidal curve; if this is not met, subject is rejected
for aa = 1:numel(prop_kitchen_mat(1,:)),
    if counter == 1,
        figure,
        hold on
    end
    subplot(x_plots,y_plots,counter)
    hold on
    %--- line fitting
    plot(my_intensity_ind(:,aa),prop_kitchen_mat(:,aa),'bo')
    p = nlinfit(my_intensity_ind(:,aa),prop_kitchen_mat(:,aa),f,[0 20/100 50/100 5/100],opts);
    new_x = linspace(0,1,total_images*10)';
    my_fit_x = f(p,new_x);
    num_unique = numel(unique(my_fit_x));
    if num_unique==1,
        while num_unique == 1,
            current_start_point = [0 20/100 50/100 5/100] .* rand;
            p = nlinfit(my_intensity_ind(:,aa),prop_kitchen_mat(:,aa),f,current_start_point,opts);
            my_fit_x = f(p,new_x);
            num_unique = unique(my_fit_x);
        end
    end
    plot(new_x,my_fit_x,'r');
    %--- average bathroom
    if range(my_fit_x) < min_prop && numel(unique(my_fit_x)) < (min_rooms_per_grouping*3) || my_fit_x(1) > my_fit_x(numel(my_fit_x)),
        text(.02,.5,'NOT ENOUGH')
        text(.02,.3,'MARGIN')
    else %we're good to go... split margin into three equal groupings; use a linear line to divy things up
        
        lin_fit = linspace(min(my_fit_x),max(my_fit_x),total_images);
        %--- low group
        low_ind = lin_fit(1:total_images/3);
        sub_average_bathrooms(aa) = max(low_ind);
        low_ind = my_fit_x(my_fit_x<=sub_average_bathrooms(aa));
        %--- high group
        high_ind = lin_fit(total_images/3*2+1:total_images);
        sub_average_kitchens(aa) = min(high_ind);
        high_ind = my_fit_x(my_fit_x>sub_average_kitchens(aa));
        %--- mid group
        mid_ind = lin_fit(total_images/3+1:total_images/3*2);
        mid_ind = my_fit_x(my_fit_x>max(low_ind) & my_fit_x<=max(mid_ind));
        if isempty(mid_ind),
            %--- this is a curve where there's a dramatic
            %difference between low and high...
            mid_ind = lin_fit(lin_fit>max(low_ind) & lin_fit<min(high_ind));
            sub_mids(aa) = mid_ind(round(numel(mid_ind)/2));
            plot(new_x(ismember(my_fit_x,max(low_ind))>0),sub_mids(aa),'ks--')
        else
            %--- Midpoint
            sub_mids(aa) = mid_ind(round(numel(mid_ind)/2));
            plot(new_x(ismember(my_fit_x,sub_mids(aa))>0),sub_mids(aa),'ks--')
        end
        plot(new_x(ismember(my_fit_x,low_ind)>0),low_ind,'-*','Color',[1 .2 0])
        plot(new_x(ismember(my_fit_x,high_ind)>0),high_ind,'-*','Color',[0 .5 1])
    end
    %---
    legend('off')
    xlabel('Proportion of Kitchen Spaciousness');
    ylabel('Proportion Kitchen Judgments');
    title(sprintf('Subject %i\rMidpoint = %.2f',aa,sub_mids(aa)))
    axis([-0.1 1.1 -0.1 1.1])
    counter = counter + 1;
    if counter > x_plots*y_plots,
        counter = 1;
    end
    
end


%--- Cross Subject performance
mean_kitchen_mat = mean(prop_kitchen_mat,2);
boot_kitchen_mat = nan(2,numel(mean_kitchen_mat));
if numel(prop_kitchen_mat(1,:)) > 1, %if we have > 1 subject
    for mm = 1:numel(mean_kitchen_mat),
        boot_kitchen_mat(1:2,mm)=bootci(10000,{@mean,prop_kitchen_mat(mm,:)},'type','bca');
    end
else
    mean_kitchen_mat = prop_kitchen_mat;
end
figure,
hold on
plot(intensity_ind(:,aa),mean_kitchen_mat,'k')
plot(intensity_ind(:,aa),boot_kitchen_mat(1,:),'r--'),
plot(intensity_ind(:,aa),boot_kitchen_mat(2,:),'r--'),
axis([-0.1 1.1 -0.1 1.1])
xlabel('Proportion of Kitchen Spaciousness');
ylabel('Proportion Kitchen Judgments');
title({'Average Proportion of kitchen judgments across subjects','Black = mean; red = 90% CI'})
fprintf('\r\r\r***************\r\r\rFinished Analyzing Data')

%---
switch produce_fmri_fits
    case 'on'
        fprintf('\r\r\r***************\r\r\rStoring Midpoints')
        out_dir = fullfile(home_dir,'fmri_midpoints');
        if ~exist(out_dir,'dir'),
            mkdir(out_dir)
        end
        my_trim_file = regexp(myfile,'\.','split');
        my_trim_file = my_trim_file{1};
        for s_idx = 1:numel(dbs_to_execute),
            sub_mid = sub_mids(aa);
            sub_average_kitchen = sub_average_kitchens(aa);
            sub_average_bathroom = sub_average_bathrooms(aa);
            save(fullfile(out_dir,sprintf('%s_%i.mat',my_trim_file,dbs_to_execute(s_idx))),'sub_mid',...
                'sub_average_kitchen','sub_average_bathroom')
        end
end
switch generate_dbs_for_each
    case 'on'
        fprintf('\r\r\r***************\r\r\rGenerating DBseqs for each sub')
        my_trim_file = regexp(myfile,'\.','split');
        my_trim_file = my_trim_file{1};
        for s_idx = 1:numel(dbs_to_execute), %subject loop
            for r_idx = 1:settings.num_runs, %run loop -- change block settings with settings.seshlimit
                settings.out_name = sprintf('%s_%i_%i',my_trim_file,dbs_to_execute(s_idx),r_idx);
                settings.out_dir = fullfile(home_dir,'db_sequences',sprintf('%s_%i',my_trim_file,dbs_to_execute(s_idx)));
                script_file = strcat(settings.out_name,'.sh');
                create_db_script(script_file,settings);
            end
        end
end
switch execute_dbs_for_each
    case 'on'
        fprintf('\r\r\r***************\r\r\rExecuting DB scripts for each sub')
        for s_idx = 1:numel(dbs_to_execute), %subject loop
            for r_idx = 1:settings.num_runs, %run loop -- change block settings with settings.seshlimit
                script_file = sprintf('%s_%i_%i.sh',my_trim_file,dbs_to_execute(s_idx),r_idx);
                unix(sprintf('nice -n %i sh %s &',settings.nice,fullfile(settings.script_dir,script_file)));
            end
        end
        kill_db_seq(my_trim_file,'status')
        %kill_db_seq(my_trim_file,'kill')
end

switch process_dbs
    case 'on'
        fprintf('\r\r\r***************\r\r\rProcessing DBs')
        my_trim_file = regexp(myfile,'\.','split');
        my_trim_file = my_trim_file{1};
        for s_idx = 1:numel(dbs_to_execute),
            settings.out_name = sprintf('%s_%i',my_trim_file,dbs_to_execute(s_idx));
            settings.out_dir = fullfile(home_dir,'db_sequences',settings.out_name);
            %%
            process_db_scripts(settings);
            %%
        end
end
