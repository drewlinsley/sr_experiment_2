%--- find and kill db calculations
function kill_db_seq(script_prefix,todo)
if ~exist('todo','var'),
    error('Give me a todo of ''kill'' or ''status''')
end
if sum(~strcmp(todo,'status') + ~strcmp(todo,'kill')) == 0,
    error('Give me a todo of ''kill'' or ''status''')
end
[~, results] = system('ps -ef');
results = regexp(results, '[\f\n\r]', 'split');
count = 0;
for kill_loop = 1:numel(results);
    if regexp(results{kill_loop},script_prefix),
        pid = regexp(results{kill_loop},'\s','split');
        switch todo
            case 'kill'
                a = unix(sprintf('kill %i',max(str2double(pid))));
                if a~=0,
                    error('\rCouldn''t find correct pid #\r')
                else
                    %fprintf('\rKilled your process')
                    count = count + 1;
                end
            case 'status'
                count = count + 1;
        end
    end
end

switch todo
    case 'kill'
        if count > 0,
            fprintf('\rKilled %i procs\r',count)
        else
            fprintf('\rCouldn''t find matching pid #\r')
        end
    case 'status'
        if count > 0,
            fprintf('\r%i sequences are running for %s\r',count,script_prefix)
        else
            fprintf('\rCouldn''t find matching pid #\r')
        end
        
end