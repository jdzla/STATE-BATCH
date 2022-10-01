def write_jobscript(batch_obj, idx, cwd, input_file, output_file):
    n_cpus  = batch_obj.comp_spec.get('n_cpus')
    pw_name = batch_obj.comp_spec.get('pw_name')
    mpi_command = batch_obj.comp_spec.get('mpi_command')
    if mpi_command is not None:
        command = f"{mpi_command} -np {n_cpus} ./{pw_name} < {input_file} > {output_file}"
    else:
        command = f"./{pw_name} < {input_file} > {output_file}"
    if (idx == 0):
         mode = 'w'
    else:
         mode = 'a'
    with open ("jobscript.txt", mode) as f:
        print (f"# Run system with idx {'{:04d}'.format(idx)}", file = f)
        print (f"pushd {cwd}", file = f)
        print (command, file = f)
        print ("popd", file = f)
        print (file = f)
