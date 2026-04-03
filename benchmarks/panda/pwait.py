#!/usr/bin/env python3

import sys
from time import sleep

from pandaclient import Client


jedi_task_id = int(sys.argv[1])
num_query_errors = 0

try:
    while True:
        status, result = Client.getTaskStatus(jedi_task_id)
        if status != 0:
            print(f"Job query error: {result}")
            num_query_errors += 1
            if num_query_errors >= 10:
                sys.exit(2)
        print(result)
        if result == "failed":
            sys.exit(1)
        elif result == "done":
            sys.exit(0)
        sleep(60)
except KeyboardInterrupt as e:
    print(f"Attempting to cancel job {jedi_task_id}")
    Client.killTask(jedi_task_id)
