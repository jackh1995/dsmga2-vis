import argparse
from termcolor import cprint
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-p','--problem', required=True)
parser.add_argument('-s','--start', required=False)
args = vars(parser.parse_args())

problems = {
    'mktrap': {
        'id': 1,
        'number': 10,
        'ell': 400,
        'success': 10
    },
    'ftrap': {
        'id': 2,
        'number': 10,
        'ell': 480,
        'success': 10
    },
    'cyctrap': {
        'id': 3,
        'number': 10,
        'ell': 400,
        'success': 10
    },
    'nks1': {
        'id': 4,
        'number': 100,
        'ell': 400,
        'success': 10
    },
    'spin': {
        'id': 5,
        'number': 100,
        'ell': 400,
        'success': 10
    },
    'maxsat': {
        'id': 6,
        'number': 1000,
        'ell': 100,
        'success': 10
    },
    'maxcut': {
        'id': 7,
        'number': 50,
        'ell': 100,
        'success': 10
    },
}

problem = args['problem']

if problem not in problems.keys():
    raise Exception(f'Problem \"{args["problem"]}\" is not defined')

if problem in {'mktrap', 'ftrap', 'cyctrap'}:
    command = f'./sweep {problems[problem]["ell"]} {problems[problem]["success"]} {problems[problem]["id"]}'

    if args['start']:
        for i in range(int(args['start']), int(problems[problem]['number'])):
            cprint(f'[{problem} {i}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split())
    else:
        for i in range(int(problems[problem]['number'])):
            cprint(f'[{problem} {i+1}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split())

if problem in {'spin', 'maxsat', 'maxcut'}:
    command = f'./sweep {problems[problem]["ell"]} {problems[problem]["success"]} {problems[problem]["id"]}'
    
    if args['start']:
        for i in range(int(args['start']), int(problems[problem]['number'])+1):
            cprint(f'[{problem} {i}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split() + [f'{i}'])
    else:
        for i in range(1, int(problems[problem]['number'])+1):
            cprint(f'[{problem} {i}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split() + [f'{i}'])

if problem == 'nks1':
    command = f'./sweep {problems[problem]["ell"]} {problems[problem]["success"]} {problems[problem]["id"]} 1'

    if args['start']:
        for i in range(int(args['start']), int(problems[problem]['number'])+1):
            cprint(f'[{problem} {i}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split() + [f'{i}'])
    else:
        for i in range(int(problems[problem]['number'])):
            cprint(f'[{problem} {i+1}/{problems[problem]["number"]}]', 'red' ,attrs=['bold'])
            subprocess.run(command.split() + [f'{i}'])