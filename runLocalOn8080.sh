originaldir=`pwd`
./bin/InterSpec --docroot . --http-address 127.0.0.1 --http-port 8080 -c ./data/config/wt_config_localweb.xml --accesslog=- --no-compression
