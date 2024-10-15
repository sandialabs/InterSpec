FROM quay.io/pypa/manylinux2014_x86_64:latest

# Install necessary packages
RUN yum update \
    && yum install -y npm \
    && npm install uglify-js -g \
    && npm install uglifycss -g \
    && npm install cmake-js -g \
    && npm install --save-dev node-addon-api --arch=x64 \
    && npm install node-api-headers \
    && npm install electron --arch=x64 \
    && npm install electron-packager

CMD ["bash"]