squeue -u $USER --format="%.10i %.20j %.2P %.2t %.7M %.16e %.3C %.7m %R" -S "${1:P,t,-p}"
running=$(squeue -u $USER -h -t running -r | wc -l)
pending=$(squeue -u $USER -h -t pending -r | wc -l)
echo "  $running jobs running; $pending more queued"
