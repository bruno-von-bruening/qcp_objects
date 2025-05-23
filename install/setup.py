#!/usr/bin/env
import os, yaml

home=os.path.realpath('..')

TAG='PROPERTY_OBJECTS'

def setup_env_vars():
    env={f'{TAG}_HOME':home}
    return env

def  env_setup_instructions(env):
    instructions=[
        f"export {k}=\"{v}\""
        for k, v in env.items()
    ]
    with open('add_to_shell_init.txt','w') as wr:
        wr.write('\n'+'\n'.join([f'# paths for {TAG} pip module']+instructions)+'\n\n')

    indent=4*' '
    instructions='\n'.join([f"Add the following commands to your shell initialization:"]+[ indent +x for x in instructions])

    # COuld also be done by just runnig cat and appending to right script
    print(instructions)

def drop_to_yaml(whole_dic):
    env_file='env_config.yaml'
    with open(env_file, 'w') as wr:
        yaml.safe_dump(whole_dic, wr)
        string=yaml.safe_dump(whole_dic)

def setup():

    env=setup_env_vars()
    env_setup_instructions(env)
    drop_to_yaml(env)

if __name__=='__main__':
    setup()
