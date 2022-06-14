"""Module containing cosine calculation code

Authors: Arjan Draisma
"""

import logging
import multiprocessing

from sklearn.metrics.pairwise import cosine_similarity

def cosine_worker(
    working_q,
    output_q,
    features
):
    """Worker function for cosine distance calculation"""
    # loop infinitely
    while True:
        # get a task
        bgc_id_a, bgc_id_b, group = working_q.get(True)
        # if the first id is none, interpret as stop signal
        if bgc_id_a is None:
            break

        # start calculation
        feature_set = features.loc[[bgc_id_a, bgc_id_b]]
        similarity = cosine_similarity(feature_set)[0,1]

        # convert to distance
        distance = 1 - similarity

        # send results on queue
        output_q.put((bgc_id_a, bgc_id_b, group, distance))
    return

def get_corr_cosine_dists(
    run,
    pairs,
    features
):
    """Function to calculate cosine distances between pairs based on a set
    of features

    Inputs:
    - run: the run object for the current run
    - pairs: a list of tuples in the format (bgc_id_a, bgc_id_b, class)
    - features: a matrix containing features per bgc id to calculate cosine
      distances on
    """
    # calculate cosine distance
    cosine_dist_corr = []

    # get the number of comparisons
    comparisons = len(pairs)

    # assign number of threads for this operation
    # in this case we want to make sure all cores are always occupied
    num_threads = run.options.cores * 4


    # create queues
    working_q = multiprocessing.Queue(num_threads)

    output_q = multiprocessing.Queue(comparisons)

    processes = []

    # start a new process for each thread
    for thread in range(num_threads):
        thread_name = "thread_cosine_" + str(thread)
        new_process = multiprocessing.Process(
            name=thread_name,
            target=cosine_worker,
            args=(
                working_q,
                output_q,
                features
            ))
        processes.append(new_process)
        new_process.start()

    # run while there are tasks
    index = 0
    done = 0
    while True:
        # set important conditions
        all_tasks_put = index == comparisons
        all_tasks_done = done == comparisons

        # break if nothing needs to be done and all results are in
        if all_tasks_put and all_tasks_done:
            break

        # if there is a spot in the queue, place a new task and increment idx
        if not working_q.full() and not all_tasks_put:
            working_q.put(pairs[index])
            index += 1
            if not working_q.full():
                continue

        # if there is a result available, retrieve and store in list
        if not output_q.empty():
            bgc_id_a, bgc_id_b, group, distance = output_q.get()
            done += 1

            cosine_dist_corr.append([bgc_id_a, bgc_id_b, group, distance])

            # print progress
            if comparisons <= 10:
                if done % int(comparisons/10) == 0:
                    log_line = "%d%% done", round(done/comparisons * 100, 1)
                    logging.debug(log_line)

    # once we are out of the loop,
    for thread_num in range(num_threads * 2):
        logging.debug("Stopping thread (%d/%d)", thread_num+1, num_threads)
        working_q.put((None, None, None))

    for process in processes:
        process.join()
        thread_name = process.name

    return cosine_dist_corr
