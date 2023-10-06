FROM mas.maap-project.org/root/maap-workspaces/base_images/vanilla:v3.0.1
ENV MAAP_PGT=PGT-2040-r0PzQ3hHeTjPkdR-vjO9qeHZIbfLJejo5XAgUsT0NegKnrvOZ1cTgQpRwDq2ylTNA1o-3604f64e67ce
COPY ./maap_runtime/submit-dps-job.py submit-dps-job.py

# ENTRYPOINT arguments get dynamically passed as part of the `docker run <name> <args>` command
# example: docker run test hello-world-aimee_ubuntu main aimeeb maap-dps-worker-16gb --params '{"one":"two"}'
ENTRYPOINT ["python", "submit-dps-job.py"]
