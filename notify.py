#!/usr/bin/env python3
import sys
import os
import smtplib
from email.message import EmailMessage
from email.utils import make_msgid, formatdate
from datetime import datetime

RECIPIENT = "mcramos22@up.edu.ph"
SENDER = "vasp-cluster@cluster.pssclabs.com"
GMAIL_MX = "142.251.10.27"  # aspmx.l.google.com for up.edu.ph

def send_notification(job_name, job_id, exit_code, nodelist):
    status_str = "SUCCESS" if str(exit_code) == "0" else f"FAILED/STOPPED (code {exit_code})"
    subject = f"[VASP {status_str}] {job_name} (Job {job_id}) on {nodelist}"
    
    # Check for last lines of OSZICAR / OUTCAR in current working directory
    last_oszicar = ""
    if os.path.exists("OSZICAR"):
        try:
            with open("OSZICAR", "r", errors="ignore") as f:
                lines = f.readlines()
                last_oszicar = "".join(lines[-10:])
        except Exception:
            pass

    body = f"""Hello Mark,

Your VASP calculation on the cluster has finished.

Job Name   : {job_name}
Job ID     : {job_id}
Status     : {status_str}
Node       : {nodelist}
Directory  : {os.getcwd()}
Timestamp  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

--- Last OSZICAR output ---
{last_oszicar if last_oszicar else 'No OSZICAR entries available.'}

Best regards,
VASP Compute Cluster
"""

    msg = EmailMessage()
    msg.set_content(body)
    msg["Subject"] = subject
    msg["From"] = SENDER
    msg["To"] = RECIPIENT
    msg["Date"] = formatdate(localtime=True)
    msg["Message-ID"] = make_msgid(domain="cluster.pssclabs.com")

    try:
        server = smtplib.SMTP(GMAIL_MX, 25, timeout=15)
        server.ehlo("cluster.pssclabs.com")
        server.send_message(msg)
        server.quit()
        print(f"[notify.py] Email notification sent to {RECIPIENT}")
    except Exception as e:
        print(f"[notify.py] Warning: Failed to send email: {e}", file=sys.stderr)

if __name__ == "__main__":
    job_name = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("SLURM_JOB_NAME", "Unknown")
    job_id = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("SLURM_JOB_ID", "Unknown")
    exit_code = sys.argv[3] if len(sys.argv) > 3 else "0"
    nodelist = sys.argv[4] if len(sys.argv) > 4 else os.environ.get("SLURM_NODELIST", "n003")
    send_notification(job_name, job_id, exit_code, nodelist)
