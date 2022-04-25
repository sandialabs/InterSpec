FROM alpine:latest

WORKDIR /mnt/host

SHELL ["/bin/bash", "-c"]

EXPOSE 8078

# Then numeric group/user value of 280 was chosen randomly; it doesnt conflict with existing groups/users on dev or public server, and is below 1000 (e.g., a system user without a home directory or default shell)
#RUN groupadd --gid 280 interspec && useradd --uid 280 --gid interspec interspec
#RUN addgroup -S interspec && adduser --disabled-password --no-create-home -S interspec -G interspec
#USER interspec
# Or just use user guest
#USER guest

# Copy our app into the container
# --chown=fullspec:fullspec
COPY build_alpine/interspec_install /var/opt/interspec

WORKDIR /var/opt/interspec

ENTRYPOINT ["/var/opt/interspec/InterSpec", "--docroot=/var/opt/interspec/html_root/", "--http-address=0.0.0.0", "--http-port=8078", "--config=/var/opt/interspec/html_root/data/config/wt_config_localweb.xml", "--accesslog=-", "--no-compression", "--userdb /var/tmp/user_data.sqlite3", "--static-data-dir", "/var/opt/interspec/html_root/data/"]


# If build fails, issue this command
# docker rmi $(docker images -f “dangling=true” -q)
# docker system prune