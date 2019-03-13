from .const import here
from subprocess import Popen, PIPE

def capture_stdouterr(command, cwd = here):
    p = Popen(command, cwd=cwd, stdout=PIPE, stderr=PIPE).communicate()
    p_out = p[0].decode('utf-8').strip()
    p_err = p[1].decode('utf-8').strip()
    return (p_out, p_err)

