
docker build -t interspec_dev_env:1.0 -f "setup_build_env.dockerfile" ../..

docker build -t makeinterspec:v1 -f "build_interspec.dockerfile" ../..

docker run  --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -i -t -v /Users/wcjohns/rad_ana/InterSpec:/interspec makeinterspec:v1 sh

#to debug: docker run  --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -i -t -v /Users/wcjohns/rad_ana/InterSpec:/interspec makeinterspec:v1 sh