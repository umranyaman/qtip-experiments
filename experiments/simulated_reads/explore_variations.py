"""

"""

import os
import sys
import logging
import shutil
from marcc_out import write_slurm

join = os.path.join
mem_gb = 8
hours = 8


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def handle_dir(dr, start_from, global_name, base_args, exp_names, exp_qsim_args, targets, submit_fh,
               use_scavenger=False):
    """
    Maybe this just creates a whole series of new Makefiles with only the SUBSAMPLING_ARGS line different?
    Then maybe the
    """
    for name, ar in zip(exp_names, exp_qsim_args):
        nm = '.'.join([global_name, name])
        new_makefile_base = '.'.join(['Makefile', global_name, name])
        logging.info('  Creating new Makefile: %s' % join(dr, new_makefile_base))
        with open(join(dr, new_makefile_base), 'w') as mk_out:
            for ln in open(join(dr, 'Makefile')):
                # 2 things to do: change the args passed to qsim and change the .out target names
                if ln.startswith('SUBSAMPLING_ARGS'):
                    mk_out.write('SUBSAMPLING_ARGS=%s %s\n' % (' '.join(base_args), ' '.join(ar)))
                else:
                    mk_out.write(ln.replace('.out', '.%s.out' % nm).replace(',out', ',%s.out' % nm))
        for fulltarget in targets:
            targdir, rule = fulltarget.split('/')
            if targdir != dr:
                continue
            orig_rule = rule
            rule = rule.replace('.out', '.%s.out' % nm)
            logging.info('    Adding job to make target: %s/%s' % (dr, rule))
            if start_from == 'inputalign':
                dest_dir = join(dr, rule)
                src_dir = join(dr, orig_rule)
                mkdir_quiet(dest_dir)
                assert os.path.exists(src_dir)
                assert os.path.exists(join(src_dir, 'input.sam'))
                logging.info('      Copying %s to new target dir' % (join(src_dir, 'input.sam')))
                shutil.copy(join(src_dir, 'input.sam'), dest_dir)
                assert os.path.exists(join(dest_dir, 'input.sam'))
            fn = '.' + rule + '.sh'
            write_slurm(rule, fn, dr, mem_gb, hours,
                        ncores=1 if start_from == 'inputalign' else 8,
                        makefile=new_makefile_base,
                        use_scavenger=use_scavenger)
            submit_fh.write('pushd %s && sbatch %s && popd\n' % (dr, fn))


def go(args, global_qsim_args, exp_names, exp_qsim_args, targets):
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S', level=logging.DEBUG)

    logging.info('Global qsim args: ' + str(global_qsim_args))
    logging.info('Experiments: ' + str(list(zip(exp_names, exp_qsim_args))))
    logging.info('Targets: ' + str(targets))

    target_dirs = set()
    for fulltarget in targets:
        _dir, _ = fulltarget.split('/')
        target_dirs.add(_dir)

    logging.info('Target dirs: ' + str(list(target_dirs)))

    with open(args.name + '_submit.sh', 'w') as submit_fh:
        # Descend into subdirectories looking for Makefiles
        for dirname, dirs, files in os.walk('.'):
            dirname = dirname.split('/')[-1]
            assert not (dirname in target_dirs and 'IGNORE' in files)
            if 'Makefile' in files and 'IGNORE' not in files and dirname in target_dirs:
                logging.info('Found a relevant Makefile: %s' % join(dirname, 'Makefile'))
                handle_dir(dirname, args.start_from, args.name, global_qsim_args, exp_names,
                           exp_qsim_args, targets, submit_fh, use_scavenger=args.use_scavenger)


def add_args(parser):
    parser.add_argument('--name', metavar='str', type=str, required=True,
                        help='Name for this overall series of experiments')
    parser.add_argument('--start-from', metavar='str', type=str, default="beginning",
                        help='What step to resume after.  options: "beginning" or "inputalign"')
    parser.add_argument('--use-scavenger', action='store_const', const=True, default=False,
                        help='Use the MARCC scavenger queue')


def parse_qsim_parameters_from_argv(argv):
    """
    :param argv: command line parameters
    :return: tuple, with 1st elt being list of parameters to this script and
    second being list of lists, each element being parameters for one qsim run
    """
    argv = argv[:]
    sections = [[], [], [], [[]], []]
    nested = [False, False, False, True, False]
    section_i = 0
    # maybe move these params into an input file
    for arg in argv:
        if arg == '==':
            section_i += 1
            continue
        if nested[section_i]:
            if arg == '--':
                sections[section_i].append([])
            else:
                sections[section_i][-1].append(arg)
        else:
            sections[section_i].append(arg)
    return tuple(sections)


if __name__ == "__main__":

    import argparse

    _parser = argparse.ArgumentParser(description='')

    add_args(_parser)

    # Some basic flags
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')

    _argv, _global_qsim_args, _exp_names, _exp_qsim_args, _targets = parse_qsim_parameters_from_argv(sys.argv)
    _args = _parser.parse_args(_argv[1:])

    go(_args, _global_qsim_args, _exp_names, _exp_qsim_args, _targets)