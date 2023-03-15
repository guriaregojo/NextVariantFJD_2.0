# Buid the image on top of another base image 

FROM ubuntu



# Commands to be executed 

# Install bcftools
RUN apt update && apt upgrade && apt-get install -y wget 


RUN wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /bin/bs

RUN chmod +x /bin/bs

RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

# Add  the necessary files 

CMD /bin/bash

#docker build -t basespace -f Dockerfile .


#docker run -it -e "HOME=${HOME}" -v $HOME:$HOME -u $(id -u):$(id -g) basespace bs auth
#docker run -it -v $HOME:/root -u $(id -u):$(id -g) basespace bs list project 

 