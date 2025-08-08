The simplest way to install is to use the ```qcp_installer``` package
(https://github.com/bruno-von-bruening/qcp_installer)

After installing the corresponding conda environment to ```qcp_installer```.
```
conda activate qcp_installer
# Go to install/ directory in this package
qcp_setup_env.py installer.yaml
```

Add the path of this repo to your local machine path.
I would recommend to generate a $HOME/.qcp_paths.sh file with all the paths and sourcing this file in the rc file of
your shell.
```
# go to the root of the git repo
echo "export QCP_OBJECTS_HOME=$(pwd)" >> $HOME/.qcp_paths.sh
```
