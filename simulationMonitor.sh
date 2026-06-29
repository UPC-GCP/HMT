#!/bin/bash

#CONTROL
MFILE="simulation_metrics.txt"
LFILE="simulation_log.txt"
counter=0
START_TIME=$(date +%s)

# HEADER
rm -f "$MFILE"
echo "MONITOR STARTED ... $(date)" > "$MFILE"

# LOOP
while pgrep -f "NSSolver" > /dev/null; do

	# METRICS
	echo "CONTROL: $(date '+%H:%M:%S')" >> "$MFILE"
	free -h | grep "Mem:" | awk '{print "RAM Used: " $3 " / " $2}' >> "$MFILE" # RAM
	uptime | awk -F'load average:' '{print "Load Avg:" $2}' >> "$MFILE" # CPU
	lscpu | grep "CPU MHz" | head -n 1 >> "$MFILE" # CLOCK
	df -h / | awk 'NR==2 {print "Disk Used: " $4 " / " $2 " (" $5 " free)"}' >> "$MFILE" # DISK
	echo "" >> "$MFILE"
	
	# TIME
	CURRENT_TIME=$(date +%s)
	DIFF_TIME=$((CURRENT_TIME - START_TIME))
	RUNTIME=$(printf '%02dh:%02dm:%02ds' $((DIFF_TIME/3600)) $((DIFF_TIME%3600/60)) $((DIFF_TIME%60)))

	# PUSH
	if [[ -n $(git status --porcelain "$LFILE") ]]; then
		git add "$MFILE" "$LFILE"
		git commit -m "Automatic: Log Update | Runtime: $RUNTIME" -q
		git push origin main -q
	else
		git add "$MFILE"
		git commit -m "Automatic: Metrics | Runtime: $RUNTIME" -q
		git push origin main -q
	fi
	
	# CONTROL
	sleep 900

done

# TIME
CURRENT_TIME=$(date +%s)
DIFF_TIME=$((CURRENT_TIME - START_TIME))
RUNTIME=$(printf '%02dh:%02dm:%02ds' $((DIFF_TIME/3600)) $((DIFF_TIME%3600/60)) $((DIFF_TIME%60)))

# FINAL PUSH
echo "MONITOR ENDED ... $(date)" >> "$MFILE"
echo "TOTAL RUNTIME: $RUNTIME" >> "$MFILE"


git add "$MFILE" "$LFILE"
git commit -m "Final Sync: Process ended." -q
git push origin main -q

# COOLDOWN
echo "Cooldown period ..."
sleep 120

# SLEEP
echo "Entering sleep mode ..."
rundll32.exe powrprof.dll,SetSuspendState 0,1,0

