
def find_path_edges(branch, tree_edges):
    branch_edges = []
    for edge in tree_edges:
        if (edge[0] in branch) and (edge[1] in branch):
            branch_edges.append(edge)
    return set(branch_edges)


def get_tree_edges(tree_paths):
    all_edges = list()
    for path in tree_paths:
        if len(path) == 2:
            all_edges.append(tuple(path))
        else:
            for i in range(len(path) - 1):
                all_edges.append((path[i], path[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def get_length_from_name(segment):
    e = int(segment.split("_")[-1])
    s = int(segment.split("_")[-2])
    return e - s


def print_logo():
    print(
    """#
    _____ __    _____ _____ _____ _____
    |  _  |  |  |  _  |  _  |     |  _  |
    |     |  |__|   __|     |   --|     |
    |__|__|_____|__|  |__|__|_____|__|__|
    /\\⌒⌒⌒/\\
    (⦿   ⦿)
    ( 'Y' )
     (   )
     (   )
     (   )
     (~ ~~~~~~~~~~)
     ( ~ ~~   ~~  )
     ( ~  ~ ~  ~  )
     (~  ~~~~~   ~)
     │ │      │ │
     │ │      │ │
    """)
