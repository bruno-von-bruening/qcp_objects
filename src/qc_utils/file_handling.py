from . import *

# Pydantic classes for file and dictionary!
def check_dir(string):
    assert os.path.isdir(string),f"Not a valid directory: {string}"
    return string
def check_file(string):
    if  os.path.isfile(string):
        return string
    elif os.path.islink(string):
        return string
    else:
        raise ValueError(f"Not a valid file: {string}")
directory=Annotated[str, BeforeValidator(check_dir)]
file=Annotated[str, BeforeValidator(check_file)]

@validate_call
def run_shell_command(
    cmd:str, 
    stdout_file: str|None='shell_output.out' , stderr_file: str|None='shell_output.err', 
    uncritical:bool=False, time_limit: float|None=None
) -> Tuple[str, str]:
    """
        time_limit: time in seconds after which to abort execution automatically (not implemented yer)
    """

    com=f"There are already redirections (\'>\') in your command, provide them with designated keywords"
    if not stdout_file is None or stderr_file is not None:
        assert not '>' in cmd, com

    if stdout_file is not None:
        cmd=cmd+f" 1> {stdout_file}"
    if stderr_file is not None:
        cmd=cmd+f" 2> {stderr_file}"
    @validate_call
    def get_out(base_str:str, the_file:file) -> str:
        with open(the_file,'r') as rd:
            content=rd.read()
        content=' | '.join([base_str, content])
        return content
        
    # Execute
    p=sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    if time_limit is None:
        stdout, stderr = [x.decode('utf-8') for x in p.communicate() ]
        if stdout_file is not None:
            stdout=get_out(stdout, stdout_file)
        if stderr_file is not None:
            stderr=get_out(stderr, stderr_file)
    else:
        raise Exception(f"Time limit not implemented yet")

    # Process the return
    if p.returncode!=0:
        if not uncritical:
            raise Exception(f"The command \'{cmd}\' failed:\n{stderr}")
    return stdout, stderr


def link_file(source:file, target:str|None=None) -> file:
    if target is None:
        target='ln.'+os.path.basename(source)
    cmd=f"ln -s {source} {target}"
    run_shell_command(cmd)
    return target

# Temporary file
class choice(str, Enum):
    tmp = 'temporary'
    src = 'source'
class temporary_file_base(BaseModel):
    """ Source is the original file, tmp is a copy of it and which specifies which file to use
    """
    source: file
    tmp_file: file | None
    which:  choice = choice.tmp
    @property
    def file(self):
        if self.which==choice.tmp:
            assert self.tmp_file is not None, f"Excepected to use temporary file but none has been define"
            return self.tmp_file
        elif self.which==choice.src:
            return self.source
    def check(self):
        assert os.path.isfile(self.source), f"Source file corrupted ({self.source} is not a file): quiting for safety"
        assert os.path.realpath(self.source)!=os.path.realpath(self.tmp_file), f"error in this class tmp_file is source file, dangerous since source could be destroyed"
        return True

class temporary_file(temporary_file_base):
    """
    Recognizes if file is linked or compressed file. 
    If not:
        works with original file
    Else:
        works with uncomprsessed copy of the file
    """
    def __init__(self, fi):
        # In case its a symlink copy to local dir
        source=fi
        tmp_file=None
        info=self.handle_source(source)
        super().__init__(**info)

    def handle_source(self,source):
        # If its a link unlink first
        tmp_file=source
        if os.path.islink(tmp_file):
            old_file=tmp_file
            tmp_file='tmp.'+os.path.basename(tmp_file)
            run_shell_command(f"cp {os.path.realpath(old_file)} {tmp_file}")
            assert os.path.isfile(tmp_file), f"Expected but does not exist: {tmp_file} (pwd: {os.getcwd()})"
            which=choice.tmp            

        ## Unpack 
        if tmp_file.endswith('.gz'):
            fi_zipped=tmp_file
            p=sp.Popen(f"gunzip {fi_zipped}", shell=True)
            p.communicate()
            tmp_file=fi_zipped.replace('.gz','')
            assert os.path.isfile(tmp_file), f"Not a file: {tmp_file}"
        if tmp_file.endswith('.xz'):
            fi_zipped=tmp_file
            p=sp.Popen(f"xz -d {fi_zipped}", shell=True)
            p.communicate()
            tmp_file=fi_zipped.replace('.xz','')
            assert os.path.isfile(tmp_file), f"Not a file: {tmp_file}"

        # If the tmp file is the same as source, then there is no tmp_file
        if os.path.realpath(tmp_file)==os.path.realpath(source):
            tmp_file=None
            which=choice.src
        # Assert the file lies in this directory and it exists
        if tmp_file is not None:
            assert os.path.relpath(os.path.dirname(os.path.realpath(tmp_file)),os.getcwd())=='.', f"Tmp_file should be dropped in this directory, tmp_file={tmp_file}"
            assert os.path.isfile(tmp_file)

        return dict(
            source=source,
            tmp_file=tmp_file,
            which=which
        )

    def remove_tmp(self, verbose=True):
        self.check()
        if self.tmp_file is not None:
            cmd=f"rm {self.tmp_file}"
            run_shell_command(cmd)
            if verbose: print(f"Deleted temporary file {self.tmp_file}")

def get_conda_base() -> directory:
    """ Return base bath of conda setup (through conda executable)"""

    # Get the conda executable
    conda_exe_var='CONDA_EXE'
    if not conda_exe_var in os.environ.keys():
        raise Exception(f"NO environment variable \'{conda_exe_var}\'")
    conda_exe_path=os.environ[conda_exe_var]

    # Truncate the end to get base path
    assert conda_exe_path.endswith('bin/conda'), conda_exe_path
    conda_exe_path=conda_exe_path.replace('bin/conda','')

    assert os.path.isdir(conda_exe_path), f"Conda base environment assumed to be \'{conda_exe_path}\' but is not a directory!"
    return conda_exe_path


@validate_call
def get_python_from_conda_env(env_name: str, base:directory=None
    ) -> file:
    if isinstance(base, type(None)):
        base=get_conda_base()
    assert isinstance(base, str)
    assert isinstance(env_name, str)

    python=os.path.join(base, 'envs', env_name, 'bin','python')
    assert os.path.isfile(python), f"Cannot find python executable under {python}, does the environment {env_name} exist on this machine?"
    return python

@validate_call
def load_json_or_yaml(the_file:file) -> Tuple[dict, list]:
    """ """ 
    def read_the_file(the_file):
        # Read the file
        with open(the_file, 'r') as rd:
            content=rd.read()
        
        # Process the content to object
        if the_file.upper().endswith( '.JSON' ):
            content=json.loads(content)
        elif the_file.upper().endswith( '.YAML' ):
            content=yaml.save_load(content)
        else:
            raise Exception(f"Unkown ending of file (\'{the_file}\')")
        return content
    
    try:
        content=read_the_file(the_file)
    except Exception as ex:
        raise Exception(f"Could not process file (\'{the_file}\') as yaml or json file:\n error={ex}")

    return content