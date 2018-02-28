FROM openmx/openmx:master

RUN apt-get update && apt-get install -y --no-install-suggests \
			      redis-tools

RUN mkdir /rsetup
ADD ./rsetup.R /rsetup/rsetup.R
RUN Rscript /rsetup/rsetup.R

ADD ./DynamicHeritability2.RData /root/DynamicHeritability2.RData
ADD ./Rprofile_worker /root/.Rprofile
ADD concatdata.R submodels.R voxAC.R /root/

ENV OMP_NUM_THREADS 1

CMD R --no-save
