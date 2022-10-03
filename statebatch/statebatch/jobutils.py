import json
import sys
import os

def write_jobscript(batch_obj, jobinfo, jobopt=None, ext="job"):
    """write_jobscript

    Parameters
    ----------
    batch_obj : Batch class object,
        The `Batch` class object
    jobinfo : dict
        Dictionary containing {idx, cwd, input_file, output_file}
    jobopt : str,
        Selection of what kind of jobscript will be produced, [SERIAL, PARALLEL]
    ext : str, optional
        Extension string for the output script, by default "job"
    """

    # Read system information
    n_cpus  = batch_obj.comp_spec.get('n_cpus')
    pw_name = batch_obj.comp_spec.get('pw_name')
    mpi_command = batch_obj.comp_spec.get('mpi_command')
    prefix = batch_obj.comp_spec.get('prefix')
    if mpi_command is not None:
        exe_command =  f"{mpi_command} -np {n_cpus} ./{pw_name}"
    else:
        exe_command = f"./{pw_name}"

    # Command maker
    def command(exe_command, input_file, output_file):
        return f"{exe_command} < {input_file} > {output_file}"

    # Serial job script
    def write_job_serial():
        """Create serial job script"""
        with open(f"{prefix}.{ext}", 'w') as fo:
            fo.write("# SERIAL FORMAT"+'\n\n')
            for idx, val in jobinfo.items():
                fo.write(f"# Run system with idx {idx:04d}" + '\n')
                fo.write(f"pushd {val['cwd']}"+'\n')
                fo.write(command(exe_command, val['input_file'], val['output_file']) + '\n')
                fo.write("popd" + '\n\n')

    # Parallel job script
    def write_job_parallel():
        """Create parallel job script"""
        with open(f"{prefix}.{ext}", 'w') as fo:
            content = {}
            for idx, val in jobinfo.items():
                content[idx+1] = {'dir':val['cwd'], 
                                'command':command(exe_command, val['input_file'], val['output_file'])}
            json_content = json.dumps(content, indent=2)
            fo.write(json_content)

    # MAIN PROTOCOL
    if jobopt.lower() == 'serial':
        write_job_serial()
    if jobopt.lower() == 'parallel':
        write_job_parallel()