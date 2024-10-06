FROM mas.maap-project.org/root/maap-workspaces/custom_images/maap_base:v4.1.0
ARG BUILD_MAPP_PGT_ARG
ENV MAAP_PGT=$BUILD_MAPP_PGT_ARG
RUN pip install git+https://github.com/MAAP-Project/maap-py.git@develop
COPY ./maap_runtime/submit-dps-job.py submit-dps-job.py
ENTRYPOINT ["bash", "-c", "export MAAP_PGT=${MAAP_PGT} && . /opt/conda/etc/profile.d/conda.sh && conda activate fire_env && python submit-dps-job.py \"$@\"", "bash"]
