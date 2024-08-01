import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import os
from datetime import datetime, timezone, timedelta
from maap.maap import MAAP
maap = MAAP(maap_host='api.maap-project.org')

now_utc = datetime.now(timezone.utc)


def send_email(job_info: list):
    if not job_info:
        return

    sender_email = "fireatlatalerts@gmail.com"
    password = os.environ['EMAIL_PW']
    #receiver_emails = "gregcorradini@gmail.com, tempest.mccabe@nasa.gov, elijah.orland@nasa.gov"
    receiver_emails = "gregcorradini@gmail.com"

    subject = f"Job Failures Since {now_utc-timedelta(hours=1)}"

    body = []
    for job in job_info:
        job_id, cmd, tag = job
        body.append(f"{job_id} running command='{cmd}' with tag='{tag}' failed")
    body = "\n".join(body)

    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = receiver_emails
    message["Subject"] = subject
    message.attach(MIMEText(body, "plain"))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_emails.split(","), message.as_string())
        print("Email sent successfully")
    except Exception as e:
        print(f"Error: {e}")
    finally:
        server.quit()


def filter_jobs_last_hour(jobs, status):
    now_utc = datetime.now(timezone.utc)

    def job_within_last_hour(job_dict):
        job_value_dict = list(job_dict.values())[0]
        if job_value_dict['status'] != status:
            return False

        time_end_iso = job_value_dict['job']['job_info']['time_end']
        time_end_dt = datetime.fromisoformat(time_end_iso.replace('Z', '+00:00')).astimezone(timezone.utc)
        return abs(now_utc - time_end_dt) <= timedelta(hours=1)

    filtered_jobs = [
        [job_value_dict['payload_id'], job_value_dict['job']['params']['_command'], job_value_dict['job']['tag']]
        for job_dict in jobs
        if (job_value_dict := list(job_dict.values())[0]) and job_within_last_hour(job_dict)
    ]

    return sorted(filtered_jobs, key=lambda l: l[0])


def list_jobs():
    jobs = maap.listJobs(username='gcorradini')
    return jobs.json()


if __name__ == '__main__':
    jobs = list_jobs()
    failed_jobs = filter_jobs_last_hour(jobs['jobs'], 'job-failed')
    print(f"[ FOUND ]: {len(failed_jobs)} failed jobs")
    for job in failed_jobs:
        print(job)
    send_email(failed_jobs)

