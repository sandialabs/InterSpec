FROM quay.io/pypa/manylinux2014_x86_64:latest

# Install necessary packages
RUN yum update \
    && yum install -y npm \
    && npm install uglify-js -g \
    && npm install uglifycss -g \
    && npm install cmake-js -g

CMD ["bash"]