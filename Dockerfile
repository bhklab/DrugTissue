FROM centos:latest

RUN sh -c '/bin/echo -e "y" | yum -y update' 

RUN sh -c 'bin/echo -e "y" | yum install -y epel-release'

RUN sh -c 'bin/echo -e "y" | yum -y update'

RUN sh -c 'bin/echo -e "y" | yum install -y R git'

RUN sh -c 'bin/echo -e "y" | yum group install "Development Tools"'

RUN git clone https://github.com/bhklab/DrugTissue.git

