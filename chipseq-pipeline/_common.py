from urllib.parse import urlparse

def parse_fastq_folder(root_folder_path):
    library_run_fastqs = {}
    fastq_root_folder = pathlib.Path(root_folder_path)
    for fastq_file in sorted(fastq_root_folder.glob('**/*')):
        if not fastq_file.is_file():
            continue
        split_path = fastq_file.relative_to(fastq_root_folder).parts
        if len(split_path) != 3:
            raise Exception(
                'The fastq folder has a non-standard structure! '
                'Please, make sure that the fastq folders only contains a folder '
                'per library, each containing a folder per run, each _only_ containing '  
                'two fastq(.gz) files.')
        library, run, fastq = split_path

        library_run_fastqs.setdefault(library, {}).setdefault(run, []).append(
            str(fastq_file.absolute()))

    return library_run_fastqs


def _check_fastq_dict_structure(library_run_fastqs):
    for library, run in library_run_fastqs.items():
        if isinstance(run, dict):
            for run, fastq_files in run.items():
                if (not isinstance(fastq_files, list) or (len(fastq_files) > 2)):
                    return False
        elif isinstance(run, list):
            if (len(run) > 2):
                return False

        else:
            return False
    return True

def _parse_sra_url(sra_url):
    parsed = urlparse(sra_url)
    srr, query = parsed.path, parsed.query
    start, end = 0, -1
    if query:
        for kv_pair in query.split('&'):
            k,v = kv_pair.split('=')
            if k == 'start':
                start = v
            if k == 'end':
                end = v
    return srr, start, end

def organize_fastqs(config):
    if isinstance(config['fastq_paths'], str):
        library_run_fastqs = parse_fastq_folder(config['fastq_paths'])
    elif isinstance(config['fastq_paths'], dict):
        if not _check_fastq_dict_structure(config['fastq_paths']):
            raise Exception(
                'An unknown format for library_fastqs! Please provide it as either '
                'a path to the folder structured as "library/run/fastqs" or '
                'a dictionary specifying the project structure.')
        # fill place holder run names
        library_run_fastqs = config['fastq_paths']
        for lib in list(library_run_fastqs.keys()):
            if isinstance(library_run_fastqs[lib], list):
                library_run_fastqs[lib] = {'lane1':library_run_fastqs[k]}

        for lib in list(library_run_fastqs.keys()):
            for run in list(library_run_fastqs[lib].keys()):
                inputs = library_run_fastqs[lib][run]
                if len(inputs) == 1 and inputs[0].startswith('sra:'):
                    srr, start, end = _parse_sra_url(inputs[0])
                    library_run_fastqs[lib][run] = [
                        'downloaded_sra_fastqs/{srr}.{start}.{end}.1.fastq.gz',
                        'downloaded_sra_fastqs/{srr}.{start}.{end}.2.fastq.gz']

    else:
        raise Exception(
            'An unknown format for library_fastqs! Please provide it as either '
            'a path to the folder with the structure library/run/fastqs or '
            'a dictionary specifying the project structure.')

    return library_run_fastqs 


def _needs_downloading(fastq_files):
    if len(fastq_files) == 1 and fastq_files[0].startswith('sra:'):
        return True
    elif (fastq_files[side].startswith('sra:') 
        or fastq_files[side].startswith('http://') 
        or fastq_files[side].startswith('ftp://')):
        return True
    else:
        return False



