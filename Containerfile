FROM fedora:latest

WORKDIR /usr/local
RUN yum install -y python3 python3-pip git gcc gfortran

RUN git clone https://github.com/bkesk/mbworkbench.git mbworkbench 
RUN cd mbworkbench && pip3 install -r Requirements.txt && pip3 install -e .

CMD ["bash"]
