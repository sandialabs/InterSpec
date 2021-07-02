originaldir=`pwd`
./bin/InterSpec --docroot . --http-address 0.0.0.0 --http-port 8080 -c ./data/config/wt_config_localweb.xml --accesslog=- --no-compression
