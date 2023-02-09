FROM mas.dit.maap-project.org/root/maap-workspaces/base_images/vanilla:dit

COPY ./maap_runtime/submit-dps-job.py submit-dps-job.py

# ENTRYPOINT arguments get dynamically passed as part of the `docker run <name> <args>` command
# example: docker run test hello-world-aimee_ubuntu main aimeeb maap-dps-worker-16gb --params '{"one":"two"}'
ENTRYPOINT ["python", "submit-dps-job.py"]
