function create_db_script(script_file,settings)


if ~exist('settings','var'),
script_file = fullfile(out_dir,'test_seq.sh');
settings.in_dir = '/Users/drewlinsley/Documents/Dropbox/deBruijn_MacOSX';  
settings.trash_dir = fullfile(settings.in_dir,'db_output');
settings.out_name = 'current_sub';
settings.out_dir = fullfile(settings.in_dir,'exp_2_scripts',settings.out_name);
settings.script_dir = fullfile(out_dir,out_name);
settings.k = 6;
settings.n = 3;
settings.numbins = 6;
settings.distfile = '/Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/exp_2_distances.txt';
settings.isi = 1500;
settings.tmin = 10;
settings.tmax = 79;
settings.numiter = 50;
settings.seshlimit = 6;
end
if ~exist(settings.out_dir,'dir'),
    mkdir(settings.out_dir)
end
if ~exist(settings.trash_dir,'dir'),
    mkdir(settings.trash_dir)
end

lines = cell(36,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lines{1} = sprintf('clear');
lines{2} = sprintf('k=%i',settings.k);
lines{3} = sprintf('n=%i',settings.n);
lines{4} = sprintf('numbins=%i',settings.numbins);
lines{5} = sprintf('distfile=%s',settings.distfile);
lines{6} = sprintf('isi=%i',settings.isi);
lines{7} = sprintf('tmax=%i',settings.tmax);
lines{8} = sprintf('numiter=%i',settings.numiter);
lines{9} = sprintf('seshlimit=%i',settings.seshlimit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lines{10} = sprintf('for ((session=1; session<=seshlimit; session++))');
lines{11} = sprintf('do');
lines{12} = sprintf('	maxdetecpow=0');
lines{13} = sprintf('for (( tmin=%i; tmin<=78; tmin++))',settings.tmin);
lines{14} = sprintf('do');
lines{15} = sprintf('for (( i=1; i<=$numiter; i++ ))');
lines{16} = sprintf('do');
lines{17} = sprintf('	#echo i=$i tmin=$tmin');
lines{18} = sprintf('        %s/./debruijn -t $k $n $numbins $distfile [$tmin,$tmax] -eval $isi -debug > %s_output_$session.txt',settings.in_dir,fullfile(settings.trash_dir,settings.out_name));
lines{19} = sprintf('	detecpow=`cat %s_output_$session.txt | grep DETECTION | awk ''{ print $3 }''`',fullfile(settings.trash_dir,settings.out_name));
lines{20} = sprintf('	correlation=`cat %s_output_$session.txt | grep CORRELATION | awk ''{ print $3}''`',fullfile(settings.trash_dir,settings.out_name));
lines{21} = sprintf('        if [ $detecpow ]');
lines{22} = sprintf('        then');
lines{23} = sprintf('	        detecpow=`echo $detecpow | bc`');
lines{24} = sprintf('                compare_result=`echo "$detecpow > $maxdetecpow" | bc`');
lines{25} = sprintf('                if test $compare_result -gt 0');
lines{26} = sprintf('	        then');
lines{27} = sprintf('		        correlation=`echo $correlation | bc`');
lines{28} = sprintf('                        echo CORRELATION = $correlation');
lines{29} = sprintf('			echo DETECTION POWER = $detecpow');
lines{30} = sprintf('			maxdetecpow=$detecpow');
lines{31} = sprintf('			cp %s_output_$session.txt %s_$session.txt',fullfile(settings.trash_dir,settings.out_name),fullfile(settings.out_dir,settings.out_name));
lines{32} = sprintf('                fi');
lines{33} = sprintf('        fi');
lines{34} = sprintf('done');
lines{35} = sprintf('done');
lines{36} = sprintf('done');

script_file = fullfile(settings.script_dir,script_file);
fid = fopen(script_file,'w');
if fid == -1,
    error('Something went wrong with the fid')
end
for mm = 1:numel(lines),
    fprintf(fid,'%s\n',lines{mm});
end
fclose(fid);