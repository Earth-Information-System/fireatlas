FROM mas.maap-project.org/root/maap-workspaces/base_images/vanilla:v3.0.1
ENV MAAP_PGT=PGT-1294-yxgc33hrHOpUCLQxBO7Hh7efXUw2fYjZXpqfTxn50nZe3cGCdD2zRwu0x0zbR-bNTog-8445019b059c
COPY ./maap_runtime/submit-dps-job.py submit-dps-job.py
RUN env | grep MAAP_PGT
# ENTRYPOINT arguments get dynamically passed as part of the `docker run <name> <args>` command
# example: docker run test hello-world-aimee_ubuntu main aimeeb maap-dps-worker-16gb --params '{"one":"two"}'
ENTRYPOINT ["python", "submit-dps-job.py"]
