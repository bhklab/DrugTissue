FROM centos:latest

RUN sh -c '/bin/echo -e "y" | yum -y update' 

RUN sh -c 'bin/echo -e "y" | yum install -y epel-release'

RUN sh -c 'bin/echo -e "y" | yum -y update'

RUN sh -c 'bin/echo -e "y" | yum install -y R git'

RUN sh -c 'bin/echo -e "y" | yum group install "Development Tools"'

RUN sh -c 'bin/echo -e "y" | yum install -y emacs'

RUN yum clean all 

RUN yum -y update

RUN yum install -y libxml2-devel.x86_64 libcurl-devel openssl-devel mesa-libGL mesa-libGL-devel mesa-libGLU mesa-libGLU-devel libX11-devel

RUN yum install -y ImageMagick-c++-devel cairo-devel 

RUN yum group install -y X11

RUN yum install -y openmpi openmpi-devel

RUN yum install -y perl-CPAN

RUN git clone https://github.com/bhklab/DrugTissue.git

RUN Rscript ./DrugTissue/TEA_dependency_installer.R







