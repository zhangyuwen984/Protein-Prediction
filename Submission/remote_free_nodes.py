import time
import subprocess
import random

def get_free_node(occupied_node):
    '''Возвращает название свободной ноды, на которую возможно отправить задачу'''

#    squeue = subprocess.run('ssh shipilov.ab@calc.cod.phystech.edu "squeue"',
#               capture_output=True, shell=True, check=True, text=True).stdout
    squeue = subprocess.run(["ssh", "shipilov.ab@calc.cod.phystech.edu", "'squeue'"],
               capture_output=True, check=True, text=True).stdout
    content=list(map(lambda x: x.split(), squeue.strip().split('\n')))
    possible_nodes=list(map(lambda x: x[8], content))[1:]

    #удаляем ноды, которые ещё не выделены и у них в статусе указано "(Priority)"
    possible_nodes=list(filter(lambda element: element!='(Priority)', possible_nodes))

    number_of_processes = float('inf')

    while number_of_processes > 90:
        if set(occupied_node) == set(possible_nodes):
            occupied_node=[]
#             print('sleep & clean')
            time.sleep(30)

        node = random.choice(list(set(possible_nodes) - set(occupied_node)))
#        number_of_processes = len(subprocess.run(f'ssh shipilov.ab@calc.cod.phystech.edu "ssh {node} ps aux | grep shipilov.ab"'
#                            , capture_output=True, shell=True, check=True, text=True).stdout.strip().split('\n'))
        number_of_processes = len(subprocess.run(["ssh", "shipilov.ab@calc.cod.phystech.edu", "ssh", f"{node}", "ps", "aux", "|", "grep", "shipilo"]
                            , capture_output=True, check=True, text=True).stdout.strip().split('\n'))

        occupied_node.append(node)
    print(node, end='')
    return node

get_free_node(occupied_node=[])
