while true
do
	sudo iw wlp1s0 measurement ftm_request conf | ts '[%Y-%m-%d %H:%M:%S]'
	sleep 1
done

