while true
do
	iw wlp1s0 link | grep signal | ts '[%Y-%m-%d %H:%M:%S]'
	sleep 1
done
