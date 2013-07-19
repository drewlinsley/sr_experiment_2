clear
k=6
n=2
numbins=6
distfile=/Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/exp_2_distances.txt
isi=1500
tmax=79
numiter=50
seshlimit=6
for ((session=1; session<=seshlimit; session++))
do
	maxdetecpow=0
for (( tmin=10; tmin<=78; tmin++))
do
for (( i=1; i<=$numiter; i++ ))
do
	#echo i=$i tmin=$tmin
        /Users/drewlinsley/Documents/Dropbox/deBruijn_MacOSX/./debruijn -t $k $n $numbins $distfile [$tmin,$tmax] -eval $isi -debug > /Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/db_output/oddball_newsubs_3_17_5_output_$session.txt
	detecpow=`cat /Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/db_output/oddball_newsubs_3_17_5_output_$session.txt | grep DETECTION | awk '{ print $3 }'`
	correlation=`cat /Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/db_output/oddball_newsubs_3_17_5_output_$session.txt | grep CORRELATION | awk '{ print $3}'`
        if [ $detecpow ]
        then
	        detecpow=`echo $detecpow | bc`
                compare_result=`echo "$detecpow > $maxdetecpow" | bc`
                if test $compare_result -gt 0
	        then
		        correlation=`echo $correlation | bc`
                        echo CORRELATION = $correlation
			echo DETECTION POWER = $detecpow
			maxdetecpow=$detecpow
			cp /Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/db_output/oddball_newsubs_3_17_5_output_$session.txt /Users/drewlinsley/Documents/Dropbox/SR_grant_experiments/Experiment_2/psychophysik_norming/db_sequences/oddball_newsubs_3_17/oddball_newsubs_3_17_5_$session.txt
                fi
        fi
done
done
done
