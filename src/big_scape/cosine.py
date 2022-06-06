import logging
import multiprocessing

from sklearn.metrics.pairwise import cosine_similarity

def cosine_worker(
    working_q,
    output_q,
    features
):
    """Worker function for cosine distance calculation"""
    while True:
        bgc_id_a, bgc_id_b, group = working_q.get(True)
        if bgc_id_a is None:
            break

        # start calculation
        feature_set = features.loc[[bgc_id_a, bgc_id_b]]
        similarity = cosine_similarity(feature_set)[0,1]

        distance = 1 - similarity
        output_q.put((bgc_id_a, bgc_id_b, group, distance))
    return

def get_corr_cosine_dists(
    run,
    pairs,
    features
):
    """Function to get a """
    # calculate cosine distance
    cosine_dist_corr = []

    comparisons = len(pairs)

    num_threads = run.options.cores * 4


    working_q = multiprocessing.Queue(num_threads)

    output_q = multiprocessing.Queue(comparisons)

    processes = []

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

    index = 0
    done = 0
    while True:
        all_tasks_put = index == comparisons
        all_tasks_done = done == comparisons

        if all_tasks_put and all_tasks_done:
            break

        if not working_q.full() and not all_tasks_put:
            working_q.put(pairs[index])
            index += 1
            if not working_q.full():
                continue

        if not output_q.empty():
            bgc_id_a, bgc_id_b, group, distance = output_q.get()
            done += 1

            cosine_dist_corr.append([bgc_id_a, bgc_id_b, group, distance])

            if comparisons > 10:
                if done % int(comparisons/10) == 0:
                    if comparisons > 0:
                        logging.debug("%d%% done", round(done/comparisons * 100, 1))

    for thread_num in range(num_threads * 2):
        working_q.put((None, None, None))

    for process in processes:
        process.join()
        thread_name = process.name

    return cosine_dist_corr
